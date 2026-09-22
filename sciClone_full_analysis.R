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

# This is a batch analysis script

eligible_patients <- read_tsv("./sciclone_patient_input.tsv")
mutation_data <- read_tsv("./sciclone_mutation_input.tsv")

head(eligible_patients)
head(mutation_data)

# Verify that the mutation data contains information from every eligible patient
sample_check <- mutation_data %>%
  filter(`Patient ID` %in% eligible_patients$`Patient ID`) %>%
  group_by(`Patient ID`, `Sample Type`) %>%
  summarise(
    n_mutations = n(),
    .groups = "drop"
  )

sample_check

# Create output directories
# Store objects, tables, and plots in their own folders within ./sciclone_results/
dir.create("sciclone_results", showWarnings = FALSE)
dir.create("sciclone_results/objects", showWarnings = FALSE)
dir.create("sciclone_results/tables", showWarnings = FALSE)
dir.create("sciclone_results/plots", showWarnings = FALSE)

# Function definition for creating datasets that fit sciclone restrictions
# Parameters: 
# patient_id -> variable representing one eligible patient with primary and metastatic tumors
# mutation_data -> matrix of information regarding patients
make_sciclone_vafs <- function(patient_id, mutation_data) {
  
  patient_data <- mutation_data %>%
    filter(`Patient ID` == patient_id)
  
  primary <- patient_data %>%
    filter(`Sample Type` == "Primary") %>%
    transmute(
      chr = Chromosome,
      st = Start_Position,
      ref = t_ref_count,
      var = t_alt_count,
      vaf = VAF * 100,
      depth = t_ref_count + t_alt_count
    )
  
  metastasis <- patient_data %>%
    filter(`Sample Type` == "Metastasis") %>%
    transmute(
      chr = Chromosome,
      st = Start_Position,
      ref = t_ref_count,
      var = t_alt_count,
      vaf = VAF * 100,
      depth = t_ref_count + t_alt_count
    )
  
  list(
    primary = primary,
    metastasis = metastasis
  )
}

# Function definition for retrieving sample names
# Parameters: 
# patient_id -> variable representing one eligible patient with primary and metastatic tumors
# mutation_data -> matrix of information regarding patients
get_sample_names <- function(patient_id, mutation_data){
  mutation_data %>%
    filter(`Patient ID`== patient_id) %>%
    distinct(`Sample Type`, Tumor_Sample_Barcode) %>%
    arrange(
      factor(
        `Sample Type`,
        levels = c("Primary", "Metastasis")
      )
    ) %>%
    pull(Tumor_Sample_Barcode)
}

# Function definition for adjusting the maximum clusters variable based on number of shared mutations 
# Parameters: 
# vafs -> matrix of VAF data
# minimumDepth -> The amount of reads required, standard = 100
# requested_max -> The default amount of maximum clusters for this analysis 
get_max_clusters <- function(vafs, minimumDepth = 100, requested_max = 5) {
  usable_sites <- Reduce(
    intersect,
    lapply(vafs, function(x) {
      paste(x$chr, x$st, sep = ":")[x$depth >= minimumDepth]
    })
  )
  n_usable <- length(usable_sites)
  
  if (n_usable < 3) {
    return(NA_integer_)
  }
  return(min(requested_max, n_usable - 1))
}

test_vafs <- make_sciclone_vafs("EV-001", mutation_data)
test_vafs

# Time for iteration 
results <- list()
run_summary <- list()

for(patient_id in eligible_patients$`Patient ID`){
  cat("Running SciClone for patient:", patient_id, "\n")
  
  # Retrieve and convert data
  vafs <- make_sciclone_vafs(patient_id, mutation_data)
  vafs <- lapply(vafs, as.data.frame)

  # Get sample names
  sample_names <- get_sample_names(patient_id, mutation_data)
  
  # Determine number of usable mutations
  usable_sites <- Reduce(
    intersect, 
    lapply(vafs, function(x) {
      paste(x$chr, x$st, sep=':')[x$depth >= 100]
    })
  )
  n_usable_sites <-length(usable_sites)
  
  max_clusters <- get_max_clusters(
    vafs = vafs, 
    minimumDepth = 100,
    requested_max = 5
  )
  
  # Remove patients with too few sites to generate enough clusters
  if(is.na(max_clusters)) {
    cat("SKIPPED:", patient_id, "\n")
    
    run_summary[[patient_id]] <- data.frame(
      Patient_ID = patient_id, 
      Status = "SKIPPED",
      usable_sites = n_usable_sites, 
      maximum_clusters = NA, 
      Error = "Fewer than three usable sites"
    )
  }
  
  # Try SciClone
  result <- tryCatch(
    {
      sciClone_result <- sciClone(
        vafs = vafs,
        sampleNames = sample_names, 
        minimumDepth = 100,
        maximumClusters = max_clusters
        )
    },
    error = function(e) {
      
      cat("FAILED:", patient_id, "\n")
      cat("Error:", conditionMessage(e), "\n")
      return(NULL)
      
    }
  )
  
  # Check if SciClone succeded 
  if(is.null(sciClone_result)) {
    run_summary[[patient_id]] <- data.frame(
      Patient_ID = patient_id,
      Status = "FAILED",
      Usable_sites = n_usable_sites,
      Maximum_clusters = max_clusters,
      Error = "SciClone failed; see console output"
    )
  }
  
  # Store information
  results[[patient_id]] <- result
  cat("Finished:", patient_id, "\n")
  saveRDS(
    sciClone_result,
    file = paste0(
      "sciclone_results/objects/",
      patient_id,
      "_sciclone.rds"
    )
  )
  # Record successful run
  run_summary[[patient_id]] <- data.frame(
    Patient_ID = patient_id,
    Status = "SUCCESS",
    Usable_sites = n_usable_sites,
    Maximum_clusters = max_clusters,
    Error = NA_character_
  )
}

run_summary_table <- bind_rows(run_summary)

run_summary_table
write_tsv(
  run_summary_table,
  "sciclone_results/tables/sciclone_run_summary.tsv"
)