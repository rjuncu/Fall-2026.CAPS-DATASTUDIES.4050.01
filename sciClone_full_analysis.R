# install.packages("dyplr")
# install.packages("remotes")
# install.packages("readr")
# install.packages("devtools")
# install.packages("stringr")
# install.packages("magick")
# library(devtools)
# if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install()
# BiocManager::install(c("IRanges","limma"))
# install_github("genome/bmm")
# install_github("genome/sciClone")
library(sciClone)
library(stringr)
library(magick)

patients <- read_tsv("./sciclone_patient_input.tsv")
mutation_data <- read_tsv("./sciclone_mutation_input.tsv")

head(patients)
head(mutation_data)

ev001_sciclone <- sciClone(
  vafs = ev001_vafs,
  sampleNames = ev001_sample_names,
  minimumDepth = 100,
  maximumClusters = 3
)
