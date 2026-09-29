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

# Time for iteration -> define lists to store results 
results <- list()
run_summary <- list()

# Function definition for applying sciclone clustering onto patients
# Parameters: 
# eligible_patients -> table containing patients with over 100 depth and 3 shared mutations
# mutation_data -> high dimensional data containing mutation info on tumors 
run_sciclone <- function(eligible_patients, mutation_data) {
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
  
}

# Call the main function 
run_sciclone(eligible_patients, mutation_data)

# Write output file for summary data
run_summary_table <- bind_rows(run_summary)
run_summary_table
write_tsv(
  run_summary_table,
  "sciclone_results/tables/sciclone_run_summary.tsv"
)

# Examine test case output Patient_ID EV-001
ev001 <- readRDS("./sciclone_results/objects/EV-001_sciclone.rds")
class(ev001)
summary(ev001@vafs.merged)

# Test case output format
ev001@vafs.merged %>%
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

rds_files <- list.files("./sciclone_results/objects", pattern = "_sciclone\\.rds$", full.names = TRUE)

# Function definition for scraping data from generated .rds output files
# Parameters:
# rds__file -> variable representing a single .rds file
extracted_sci_results <- function(rds_file) {
  patient_id<- str_extract(
    basename(rds_file), 
    "EV-[0-9]+"
  )
  sc <- readRDS(rds_file)
  merged <- sc@vafs.merged
  
  primary_vaf_col <- paste0(patient_id, "-P.vaf")
  metastasis_vaf_col <- paste0(patient_id, "-M.vaf")
  primary_depth_col <- paste0(patient_id, "-P.depth")
  metastasis_depth_col <- paste0(patient_id, "-M.depth")
  
  if ("cluster.prob" %in% names(merged)) {
    cluster_prob_values <- merged$cluster.prob
  } else {
    cluster_prob_values <- rep(NA_real_, nrow(merged))
  }
  
  result <- merged %>%
    transmute(
      Patient_ID = patient_id, 
      chr = chr, 
      st = st, 
      Primary_VAF = .data[[primary_vaf_col]], 
      Metastasis_VAF = .data[[metastasis_vaf_col]], 
      Primary_depth = .data[[primary_depth_col]],
      Metastasis_depth = .data[[metastasis_depth_col]], 
      adequateDepth = adequateDepth, 
      cluster = cluster, 
      cluster_prob = cluster_prob_values
    )
  return(result)
}

# Test the extraction function 
test_result <- extracted_sci_results("./sciclone_results/objects/EV-001_sciclone.rds")
head(test_result)

# Function definition for extracting .rds files
# Parameters: 
# Function will be applied into rds_files object, using file as a variable for one input
rds_check <- lapply(
  rds_files,
  function(file) {
    
    filename_patient <- str_extract(
      basename(file),
      "EV-[0-9]+"
    )
    sc <- readRDS(file)
    object_samples <- sc@sampleNames
    
    object_patient <- str_extract(
      object_samples[1],
      "EV-[0-9]+"
    )
    
    data.frame(
      File_Patient = filename_patient,
      Object_Patient = object_patient,
      Sample_1 = object_samples[1],
      Sample_2 = object_samples[2],
      Match = filename_patient == object_patient,
      stringsAsFactors = FALSE
    )
  }
)

rds_check <- bind_rows(rds_check)

rds_check
# rds_check %>% filter(!Match)
# table(rds_check$Match)
# rds_check %>%
#   count(Object_Patient) %>%
#   filter(n > 1)
# rds_check %>%
#   filter(Object_Patient %in% c("EV-056", "EV-064"))

# Patients EV-058 and EV-066 failed sciclone clustering, so they will be excluded 
# Remove and inspect their .rds files
file.remove("./sciclone_results/objects/EV-058_sciclone.rds", "./sciclone_results/objects/EV-066_sciclone.rds")
list.files("./sciclone_results/objects/")

# Prepare for ClonEvol
rds_files <- list.files("./sciclone_results/objects", pattern = "_sciclone\\.rds$", full.names = TRUE)

# Extract all files by applying the function defined above
all_sciclone_results <- lapply(rds_files, extracted_sci_results)
clonevol_input <- bind_rows(all_sciclone_results)

# Remove samples without adequate depth
#clonevol_input <- clonevol_input %>%
#  filter(adequateDepth == 1)

# Inspect the resulting data
dim(clonevol_input)
head(clonevol_input)
summary(clonevol_input)

# Check the SciClone cluster distribution
clonevol_input %>%
  count(Patient_ID, cluster, name = "n_mutations") %>%
  arrange(Patient_ID, cluster)

# QC: inspect SciClone cluster assignments and missing cluster probabilities

clonevol_input %>%
  group_by(Patient_ID) %>%
  summarise(
    n_mutations = n(),
    n_clustered = sum(!is.na(cluster) & cluster > 0),
    n_unclustered = sum(is.na(cluster) | cluster == 0),
    n_clusters = n_distinct(cluster[!is.na(cluster) & cluster > 0]),
    n_missing_cluster_prob = sum(is.na(cluster_prob)),
    .groups = "drop"
  )

clonevol_input %>%
  filter(cluster == 0) %>%
  select(
    Patient_ID,
    chr,
    st,
    Primary_VAF,
    Metastasis_VAF,
    Primary_depth,
    Metastasis_depth,
    cluster,
    cluster_prob
  )

clonevol_input %>%
  filter(is.na(cluster)) %>%
  select(
    Patient_ID,
    chr,
    st,
    Primary_VAF,
    Metastasis_VAF,
    Primary_depth,
    Metastasis_depth,
    cluster,
    cluster_prob
  )

# Save the clustered data
sciclone_clustered <- clonevol_input %>%
  filter(cluster > 0)
sciclone_clustered_NA <- clonevol_input %>%
  filter(is.na(cluster))
sciclone_clustered_0 <- clonevol_input %>%
  filter(cluster == 0)

head(sciclone_clustered)

sciclone_clustered %>%
  count(Patient_ID, cluster) %>% 
  arrange(Patient_ID, cluster)

# Create tables
write_tsv(sciclone_clustered, "./sciclone_results/tables/sciclone_clustered")
write_tsv(sciclone_clustered_NA, "./sciclone_results/tables/sciclone_NA")
write_tsv(sciclone_clustered_0, "./sciclone_results/tables/sciclone_0")
write_tsv(clonevol_input, "./sciclone_results/tables/sciclone_total")
