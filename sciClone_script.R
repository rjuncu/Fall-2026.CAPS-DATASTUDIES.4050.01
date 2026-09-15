getwd()
#install.packages("remotes")
#install.packages("readr")
library(remotes)
library(readr)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("IRanges","limma"))

install_github("genome/bmm")
install_github("genome/sciClone")

library(sciClone)

v1 = read_tsv("annotated_mutations.tsv")
