getwd()
#install.packages("remotes")
#install.packages("readr")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# BiocManager::install(c("IRanges","limma"))
# 
# install_github("genome/bmm")
# install_github("genome/sciClone")

library(sciClone)
library(remotes)
library(readr)

#Read input data and filter for samples with a coverage depth of at least 20
raw_input <- read_tsv("annotated_mutations.tsv")
dim(raw_input)
raw_input <- raw_input %>%
  filter((t_ref_count+t_alt_count) >= 20)

#Filter input data to 5 columns
#Chromosome, genomic position, reference reads, variant reads, and VAF
input_data <- raw_input %>%
  select(Chromosome, Start_Position, t_ref_count, t_alt_count, VAF)

summary(input_data)
dim(input_data)
