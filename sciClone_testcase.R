#getwd()
#install.packages("remotes")
#install.packages("readr")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# BiocManager::install(c("IRanges","limma"))
# 
# install_github("genome/bmm")
# install_github("genome/sciClone")

#REMEMBER The goal of SciClone Analysis is: 
#Do the observed VAFs between primary and metastatic tumors from the same patient 
#form clusters that could represent different populations of mutations?

#The goal of this script is to validate the SciClone analysis on an example patient

library(sciClone)
library(remotes)
library(readr)
library(dplyr)

# For more information on SciClone
# packageVersion("sciClone")
# args(sciClone)
# getAnywhere(sciClone)

# Read input data
raw_input <- read_tsv("annotated_mutations.tsv")

# There are 781 mutation records from 138 unique tumor samples
summary(raw_input$Tumor_Sample_Barcode)

# Filter for samples with a coverage depth of at least 100
depth_100_data <- raw_input %>%
  filter((t_ref_count+t_alt_count) >= 100)

# Check for duplicate genomic positions within samples -> should return 0 
depth_100_data %>%
  count(Tumor_Sample_Barcode, Chromosome, Start_Position) %>%
  filter(n > 1)

# Investigating all samples
# Mutation_ID = chromosome:start_position
# Sorting mutations into 4 categories:
# shared, primary only, metastatic only, total unique
shared_summary <- depth_100_data %>%
  mutate(
    Mutation_ID = paste(Chromosome, Start_Position, sep = ":")
  ) %>%
  group_by(`Patient ID`, Mutation_ID) %>%
  summarise(
    in_primary = any(`Sample Type` == "Primary"),
    in_metastasis = any(`Sample Type` == "Metastasis"),
    .groups = "drop"
  ) %>%
  group_by(`Patient ID`) %>%
  summarise(
    Shared_sites = sum(in_primary & in_metastasis),
    Primary_only = sum(in_primary & !in_metastasis),
    Metastasis_only = sum(!in_primary & in_metastasis),
    Total_unique_sites = n(),
    .groups = "drop"
  )
# There are more metastatic mutations than primary mutations for various reasons
# Sequence coverage, filtering, calling sensitivity, tumor behavior changes...
# Number of mutations remaining after filtering for depth => 100 is 137 instead of 138
table(depth_100_data$`Sample Type`)

# Create a table of eligible patients (63)
# Whose primary and metastatic tumors share at least 3 mutations over 100 depth
sciclone_patient_table <- shared_summary %>%
  filter(Shared_sites >= 3) %>%
  mutate(
    Primary_Sample = paste0(`Patient ID`, "-P"),
    Metastasis_Sample = paste0(`Patient ID`, "-M")
  ) %>%
  select(
    `Patient ID`,
    Primary_Sample,
    Metastasis_Sample,
    Shared_sites,
    Primary_only,
    Metastasis_only,
    Total_unique_sites
  )
sciclone_patient_table

# Investigating which sample has less than 100 depth
all_samples <- unique(raw_input$Tumor_Sample_Barcode)

filtered_samples <- unique(depth_100_data$Tumor_Sample_Barcode)

# EV-034 PATIENT ANALYSIS 
setdiff(all_samples, filtered_samples)
# primary
raw_input %>%
  filter(Tumor_Sample_Barcode == "EV-034-P") %>%
  summarise(
    total_mutations = n(),
    max_depth = max(t_ref_count + t_alt_count, na.rm = TRUE),
    mean_depth = mean(t_ref_count + t_alt_count, na.rm = TRUE),
    mutations_depth_100 = sum((t_ref_count + t_alt_count) >= 100)
  )
 # metastatic
raw_input %>%
  filter(Tumor_Sample_Barcode == "EV-034-M") %>%
  summarise(
    total_mutations = n(),
    max_depth = max(t_ref_count + t_alt_count, na.rm = TRUE),
    mean_depth = mean(t_ref_count + t_alt_count, na.rm = TRUE),
    mutations_depth_100 = sum((t_ref_count + t_alt_count) >= 100)
  )

# EV001 PATIENT ANALYSIS
# The primary tumor for patient EV-034-P has 48 depth, but the metastatic tumor has 2877 depth
# Examine the structure of a single patient paired tumors test case
ev001 <- depth_100_data %>%
  filter(Tumor_Sample_Barcode %in% c("EV-001-P", "EV-001-M"))
ev001
# primary
ev001_primary_sc <- raw_input %>%
  filter(Tumor_Sample_Barcode == "EV-001-P") %>%
  transmute(
    chr = Chromosome,
    st = Start_Position,
    ref = t_ref_count,
    var = t_alt_count,
    vaf = VAF * 100,
    depth = (t_ref_count+t_alt_count)
  )
# metastatic
ev001_metastasis_sc <- raw_input %>%
  filter(Tumor_Sample_Barcode == "EV-001-M") %>%
  transmute(
    chr = Chromosome,
    st = Start_Position,
    ref = t_ref_count,
    var = t_alt_count,
    vaf = VAF * 100,
    depth = (t_ref_count+t_alt_count)
  )
ev001_primary_sc
ev001_metastasis_sc


#The only mutation present in the metastatic tumor and not the primary tumor
#is at position 30691872

#Observe the differences in VAF between the two tumors from the same test patient
# EV-001-P       EV-001-M
# 
# 0.157          0.245
# 0.143          0.121
# 0.197          0.240
# 0.205          0.204
# 0.193          0.240
# 0.202          0.363
# 0.256          0.408

#Increased VAF = increased percentage of reads contain mutated copies 
#Relative difference in % increase -> some mutations increase in copy number faster than others


#Create a list of two tables, with the same columns
#one containing primary, the other with metastatic data 
ev001_vafs <- list(
  ev001_primary_sc,
  ev001_metastasis_sc
)

#Earlier filtered for mutations with depth of 100 due to sciclone restrictions
#this can cause issues, for example
# Mutation_ID 3:30691872
# Primary:      not detected
# Metastasis:   detected, 205×
# potentially informative metastasis-specific event
 
# Mutation_ID 5:67589598
# Primary:      detected, 280×
# Metastasis:   only 33×
# insufficient evidence to call it absent in metastasis, 
# but not enough depth for sciclone

# Same for 
# Mutation_ID 7:140453193
# Primary:      detected, 306×
# Metastasis:   only 50×

ev001_sample_names <- c(
  "EV-001-P",
  "EV-001-M"
)

ev001_vafs

#average sequencing depth for primary and metastatic
sapply(ev001_vafs, function(x) {
  mean(x$ref + x$var)
})

ev001_vafs <- lapply(ev001_vafs, as.data.frame)

# Run SciClone on EV-001 with 3 clusters as a test
ev001_sciclone <- sciClone(
  vafs = ev001_vafs,
  sampleNames = ev001_sample_names,
  minimumDepth = 100,
  maximumClusters = 3
)

#This test case created one cluster
#Among the mutations that passed the filtering criteria (13 recorded, 8 unique, 5 high quality), 
#the VAF distributions were modeled as one cluster in the two-sample analysis.
#Cluster 1 pi = 1.000 center = 0.204 or 20.4% VAF

#One outlier was detected
#EV-001-M.vaf
#2    0.4075829
#VAF of 40.76% correlates to Mutation_ID 17:7578406
#labeled an outlier instead of a second cluster

#learn more about SciClone objects
# class(ev001_sciclone)
# slotNames(ev001_sciclone)
# ev001_sciclone@sampleNames
# ev001_sciclone@purities
# ev001_sciclone@clust

# Mutation_ID 3:30691872 is row 3 of this table 
# because this mutation was only present in the metastatic tumor, sciclone assigns NA for ref, var, vaf, and depth
# this is also happening to mutations that were recorded in the primary, but did not have over 100 depth in the metastatic
head(ev001_sciclone@vafs.merged)

# Need to distinguish between marginal and full-dimensional analysis
# Marginal is primary and metastatic separately 
# # full is them together - requires joint adequate depth criteria
# SciClone:
# Which mutations have similar VAF patterns and therefore can be grouped into clusters?
# CloneEvol:
# How can those inferred clones be arranged as an evolutionary relationship across samples/time points?

ev001_results <- ev001_sciclone@vafs.merged %>%
  select(
    chr,
    st,
    `EV-001-P.vaf`,
    `EV-001-M.vaf`,
    `EV-001-P.depth`,
    `EV-001-M.depth`,
    adequateDepth,
    cluster,
    cluster.prob
  )

ev001_results

# Now that SciClone has been successfully run on a test example, time to export data
# 2 Datasets: 
# List of eligible patients and their shared mutations
# List of mutational data associated with these patients 
sciclone_mutation_table <- depth_100_data %>%
  filter(`Patient ID` %in% sciclone_patient_table$`Patient ID`) %>%
  mutate(
    Mutation_ID = paste(Chromosome, Start_Position, sep=":")
  )

# Verify that the two datasets have the same number of patients
n_distinct(sciclone_mutation_table$`Patient ID`)
nrow(sciclone_patient_table)

# Write Files! 
write_tsv(sciclone_patient_table, "sciclone_patient_input.tsv")
write_tsv(sciclone_mutation_table, "sciclone_mutation_input.tsv")
