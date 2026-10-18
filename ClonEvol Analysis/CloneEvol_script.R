install.packages("jellyfisher")
install.packages('devtools')
library(devtools)
install_github('hdng/clonevol')
install.packages('gridBase')
install.packages('gridExtra')
install.packages('ggplot2')
install.packages('igraph')
install.packages('packcircles')
install_github('hdng/trees')
install.packages("readr")

library(jellyfisher)
library(clonevol)
library(readr)
library(dplyr)
getwd()

# This script follows the tutorial at https://raw.githubusercontent.com/hdng/clonevol/master/vignettes/clonevol.pdf


# ClonEvol analysis across all patients
# Load SciClone results
sciclone_total <- read_tsv("./sciclone_total.tsv")
summary(sciclone_total)
dim(sciclone_total)
names(sciclone_total)


# Create output directories
dir.create("clonevol_results/objects", recursive = TRUE, showWarnings = FALSE)

dir.create("clonevol_results/tables", recursive = TRUE, showWarnings = FALSE)

# Function definition for preparing patient data for clonevol bootstrapping
# Parameters: patient_data -> Represents a patient's record in sciclone_total
prepare_clonevol_input <- function(patient_data) {
  
  x <- patient_data %>%
    filter(
      adequateDepth == 1,
      !is.na(cluster),
      cluster > 0,
      !is.na(Primary_VAF),
      !is.na(Metastasis_VAF)
    ) %>%
    transmute(
      chr = chr,
      st = st,
      cluster = as.integer(cluster),
      Primary = Primary_VAF,
      Metastasis = Metastasis_VAF
    ) %>%
    distinct(chr, st, .keep_all = TRUE)
  
  # ClonEvol requires consecutive positive cluster IDs.
  if (nrow(x) > 0) {
    old_clusters <- sort(unique(x$cluster))
    x$cluster <- match(x$cluster, old_clusters)
  }
  
  x
}

# Function definition for running ClonEvol for one patient
# Parameters: 
#   patient_id <- represents a unique patient
#   patient_data <- information regarding cluster connections 
#   num_boots <- number of bootstrap resampling 
#   initiation_model <- polyclonal worked, while monoclonal setting did not
run_clonevol_patient <- function(
    patient_id,
    patient_data,
    num_boots = 1000,
    initiation_model = "polyclonal") {
  
  cat("\n----------------------------------------\n")
  cat("Patient:", patient_id, "\n")
  cat("----------------------------------------\n")
  
  # Prepare patient-specific variants
  x <- prepare_clonevol_input(patient_data)
  
  n_variants <- nrow(x)
  n_clusters <- dplyr::n_distinct(x$cluster)
  
  # Return a skipped result if there is insufficient input
  if (n_variants == 0 || n_clusters < 2) {
    return(list(
      patient_id = patient_id,
      status = "SKIPPED",
      n_variants = n_variants,
      n_clusters = n_clusters,
      n_models = NA_integer_,
      error = "No eligible variants or fewer than two clusters",
      result = NULL
    ))
  }
  
  # Validate input
  if (any(!is.finite(x$Primary)) ||
      any(!is.finite(x$Metastasis)) ||
      any(x$Primary < 0 | x$Primary > 100) ||
      any(x$Metastasis < 0 | x$Metastasis > 100)) {
    return(list(
      patient_id = patient_id,
      status = "SKIPPED",
      n_variants = n_variants,
      n_clusters = n_clusters,
      n_models = NA_integer_,
      error = "Invalid or non-finite VAF values",
      result = NULL
    ))
  }
  
  # Set up sample columns and groups
  vaf_col_names <- c("Primary", "Metastasis")
  
  sample_groups <- c(
    Primary = "P",
    Metastasis = "M"
  )
  
  # Run ClonEvol; capture errors so one patient cannot
  # stop the entire batch.
  set.seed(123)
  
  y <- tryCatch(
    infer.clonal.models(
      variants = x,
      cluster.col.name = "cluster",
      vaf.col.names = vaf_col_names,
      sample.groups = sample_groups,
      cancer.initiation.model = initiation_model,
      subclonal.test = "bootstrap",
      subclonal.test.model = "non-parametric",
      num.boots = num_boots,
      founding.cluster = NULL,
      cluster.center = "mean",
      ignore.clusters = NULL,
      min.cluster.vaf = 0.01,
      sum.p = 0.05,
      alpha = 0.05
    ),
    error = function(e) e
  )
  
  # Record errors and continue to the next patient
  if (inherits(y, "error")) {
    error_message <- conditionMessage(y)
    
    cat("FAILED:", patient_id, "\n")
    cat("Reason:", error_message, "\n")
    
    return(list(
      patient_id = patient_id,
      status = "FAILED",
      n_variants = n_variants,
      n_clusters = n_clusters,
      n_models = NA_integer_,
      error = error_message,
      result = NULL
    ))
  }
  
  # Count consensus models if available
  n_models <- if (!is.null(y$num.matched.models)) {
    y$num.matched.models
  } else {
    NA_integer_
  }
  
  # Save successful result
  output_file <- file.path(
    "clonevol_results",
    "objects",
    paste0(patient_id, "_clonevol.rds")
  )
  
  saveRDS(y, output_file)
  
  cat("SUCCESS:", patient_id, "\n")
  cat("Eligible variants:", n_variants, "\n")
  cat("Clusters:", n_clusters, "\n")
  cat("Matched models:", n_models, "\n")
  
  list(
    patient_id = patient_id,
    status = "SUCCESS",
    n_variants = n_variants,
    n_clusters = n_clusters,
    n_models = n_models,
    error = NA_character_,
    result = y
  )
}

# Function definition for iterating over all patients
# Parameters: 
#   sciclone_table <- represents input from sciclone clustering
run_clonevol_all <- function(sciclone_table,
                             num_boots = 1000,
                             initiation_model = "polyclonal") {
  
  # Check that required columns exist
  required_columns <- c(
    "Patient_ID", "adequateDepth", "cluster",
    "chr", "st", "Primary_VAF", "Metastasis_VAF"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(sciclone_table)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  patient_ids <- unique(
    na.omit(as.character(sciclone_table$Patient_ID))
  )
  
  results <- list()
  summaries <- list()
  
  for (patient_id in patient_ids) {
    
    patient_data <- sciclone_table %>%
      filter(as.character(Patient_ID) == patient_id)
    
    patient_result <- run_clonevol_patient(
      patient_id = patient_id,
      patient_data = patient_data,
      num_boots = num_boots,
      initiation_model = initiation_model
    )
    
    # Keep the full result in memory only if successful
    if (!is.null(patient_result$result)) {
      results[[patient_id]] <- patient_result$result
    }
    
    # Keep a compact summary for every patient
    summaries[[patient_id]] <- data.frame(
      Patient_ID = patient_result$patient_id,
      Status = patient_result$status,
      Eligible_variants = patient_result$n_variants,
      Number_of_clusters = patient_result$n_clusters,
      Matched_models = patient_result$n_models,
      Error = patient_result$error,
      stringsAsFactors = FALSE
    )
  }
  
  list(
    results = results,
    summary = dplyr::bind_rows(summaries)
  )
}

# Run the batch
clonevol_output <- run_clonevol_all(
  sciclone_table = sciclone_total,
  num_boots = 1000,
  initiation_model = "polyclonal"
)

clonevol_results <- clonevol_output$results
clonevol_summary <- clonevol_output$summary

# Save and inspect the summary
write_tsv(clonevol_summary, "clonevol_results/tables/clonevol_run_summary.tsv")

print(clonevol_summary)

table(clonevol_summary$Status, useNA = "ifany")



