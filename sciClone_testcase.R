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

#get info on sciclone arguments
# packageVersion("sciClone")
# args(sciClone)
# getAnywhere(sciClone)

#Read input data
raw_input <- read_tsv("annotated_mutations.tsv")

# There are 781 mutation records from 138 unique tumor samples
summary(raw_input$Tumor_Sample_Barcode)

#Filter for samples with a coverage depth of at least 100
depth_100_data <- raw_input %>%
  filter((t_ref_count+t_alt_count) >= 100)

#there are more metastatic mutations than primary mutations for various reasons
#sequence coverage, filtering, calling sensitivity, tumor behavior changes...
#number of mutations remaining after filtering for depth => 100 is 137 instead of 138
table(depth_100_data$`Sample Type`)

#investigating which sample has less than 100 depth
all_samples <- unique(raw_input$Tumor_Sample_Barcode)

filtered_samples <- unique(depth_100_data$Tumor_Sample_Barcode)

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

#the primary tumor for patient EV-034-P has 48 depth, but the metastatic tumor has 2877 depth

#Examine the structure of a single patient paired tumors test case
ev001 <- depth_100_data %>%
  filter(Tumor_Sample_Barcode %in% c("EV-001-P", "EV-001-M"))

#Clean up columns and add a mutationid and depth 
#mutationid -> address using chromosome:startposition
ev001 <- ev001[, c(
  "Tumor_Sample_Barcode",
  "Chromosome",
  "Start_Position",
  "t_ref_count",
  "t_alt_count",
  "VAF"
)]

ev001 <- ev001 %>%
  mutate(Mutation_ID = paste(Chromosome, Start_Position, sep = ":"), Depth = t_ref_count + t_alt_count)

ev001

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

ev001_primary_sc <- ev001 %>%
  filter(Tumor_Sample_Barcode == "EV-001-P") %>%
  transmute(
    chr = Chromosome,
    st = Start_Position,
    ref = t_ref_count,
    var = t_alt_count,
    vaf = VAF * 100,
    depth = Depth
  )

ev001_metastasis_sc <- ev001 %>%
  filter(Tumor_Sample_Barcode == "EV-001-M") %>%
  transmute(
    chr = Chromosome,
    st = Start_Position,
    ref = t_ref_count,
    var = t_alt_count,
    vaf = VAF * 100,
    depth = Depth
  )
ev001_primary_sc
ev001_metastasis_sc

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
# SciClone is answering approximately:
#   
#   Which mutations have similar VAF patterns and therefore can be grouped into clusters?
#   
#   CloneEvol is answering a different question:
#   
#   How can those inferred clones be arranged as an evolutionary relationship across samples/time points?

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

#Now that SciClone has been successfully run on a test example, time to find all the eligible samples in the study
#investigating all samples 
#sorting their mutations into four categories
#shared, primary only, metastatic only, total unique
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
    shared_sites = sum(in_primary & in_metastasis),
    primary_only = sum(in_primary & !in_metastasis),
    metastasis_only = sum(!in_primary & in_metastasis),
    total_unique_sites = n(),
    .groups = "drop"
  )

#Understand how many shared mutations exist in the average patient
# Validate that patient EV-034 is the only patient with no shared mutations
shared_summary %>%
  arrange(shared_sites)
shared_summary %>%
  summarise(Pairs_over_3 = sum(shared_sites >= 3))


# Check for duplicate genomic positions within samples -> should return 0 
depth_100_data %>%
  count(Tumor_Sample_Barcode, Chromosome, Start_Position) %>%
  filter(n > 1)


ev001_shared <- ev001 %>%
  group_by(Mutation_ID) %>%
  filter(n_distinct(Tumor_Sample_Barcode) == 2) %>%
  ungroup()

ev001_shared
