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
install.packages("webshot2")

library(readr)
library(dplyr)
library(clonevol)
library(jellyfisher)
library(htmlwidgets)
library(webshot2)
dir.create("jellyfisher_results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("jellyfisher_results/plots", recursive = TRUE, showWarnings = FALSE)

# Set TRUE to open each plot in the RStudio Viewer as it is made
show_in_viewer <- TRUE

# Read the SciClone and mutation annotation data
sciclone_total <- read_tsv("C:/Users/rjunc/Documents/4050_project/Fall-2026.CAPS-DATASTUDIES.4050.01/SciClone Analysis/sciclone_results/tables/sciclone_total.tsv", show_col_types = FALSE)
annotations <- read_tsv("C:/Users/rjunc/Documents/4050_project/Fall-2026.CAPS-DATASTUDIES.4050.01/SciClone Analysis/annotated_mutations.tsv", show_col_types = FALSE)

# Check that the input data has the required columns
required_sciclone <- c(
  "Patient_ID", "chr", "st", "cluster",
  "adequateDepth", "Primary_VAF", "Metastasis_VAF"
)

if (!all(required_sciclone %in% names(sciclone_total))) {
  stop("The SciClone table is missing required columns.")
}

if (!all(c("Chromosome", "Start_Position") %in% names(annotations))) {
  stop("The annotation table needs Chromosome and Start_Position columns.")
}

# Function definitions to normalize chromosome and position values for joining
normalize_chr <- function(x) {
  toupper(sub("^chr", "", trimws(as.character(x)),
              ignore.case = TRUE))
}

normalize_pos <- function(x) {
  as.character(as.numeric(x))
}

# Choose the best available annotation columns
pick_column <- function(candidates, available) {
  found <- candidates[candidates %in% available]
  if (length(found)) found[1] else NA_character_
}

gene_col <- pick_column(
  c("Hugo_Symbol", "Gene_Symbol", "Gene"),
  names(annotations)
)

protein_col <- pick_column(
  c("HGVSp_Short", "HGVSp", "Protein_Change",
    "ProteinChange", "AAChange", "Amino_acids"),
  names(annotations)
)

class_col <- pick_column(
  c("Variant_Classification", "Consequence"),
  names(annotations)
)

if (is.na(gene_col)) {
  stop("No gene-symbol column found in the annotation file.")
}

message("Gene annotation column: ", gene_col)
message("Protein-change column: ", ifelse(is.na(protein_col), "not found", protein_col))

# Combine gene name with protein change where available
annotation_lookup <- annotations %>%
  transmute(
    chr_key = normalize_chr(Chromosome),
    pos_key = normalize_pos(Start_Position),
    gene = as.character(.data[[gene_col]]),
    protein = if (!is.na(protein_col)) {
      as.character(.data[[protein_col]])
    } else {
      NA_character_
    },
    variant_class = if (!is.na(class_col)) {
      as.character(.data[[class_col]])
    } else {
      NA_character_
    }
  ) %>%
  mutate(
    mutation_label = case_when(
      !is.na(gene) & gene != "" &
        !is.na(protein) & protein != "" & protein != "NA" ~
        paste(gene, protein),
      
      !is.na(gene) & gene != "" &
        !is.na(variant_class) & variant_class != "" ~
        paste0(gene, " (", variant_class, ")"),
      
      !is.na(gene) & gene != "" ~ gene,
      
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(chr_key, pos_key) %>%
  summarise(
    mutation_label = {
      z <- unique(mutation_label[
        !is.na(mutation_label) & mutation_label != ""
      ])
      if (length(z)) paste(z, collapse = " / ") else NA_character_
    },
    .groups = "drop"
  )


# Prepare mutation labels for one patient's clones
prepare_mutation_labels <- function(patient_id) {
  
  variants <- sciclone_total %>%
    filter(
      as.character(Patient_ID) == patient_id,
      adequateDepth == 1,
      !is.na(cluster),
      as.numeric(as.character(cluster)) > 0,
      !is.na(Primary_VAF),
      !is.na(Metastasis_VAF)
    ) %>%
    mutate(
      original_cluster = as.integer(as.character(cluster)),
      chr_key = normalize_chr(chr),
      pos_key = normalize_pos(st)
    ) %>%
    distinct(chr_key, pos_key, .keep_all = TRUE)
  
  if (nrow(variants) == 0) {
    stop("No eligible SciClone variants found for ", patient_id)
  }
  
  # Match the consecutive cluster numbering used by the batch runner
  variants <- variants %>%
    mutate(
      clone_id = match(
        original_cluster,
        sort(unique(original_cluster))
      )
    ) %>%
    left_join(
      annotation_lookup,
      by = c("chr_key", "pos_key")
    ) %>%
    mutate(
      mutation_label = if_else(
        is.na(mutation_label) | mutation_label == "",
        paste0("Unannotated chr", chr, ":", st),
        mutation_label
      )
    )
  
  # One row per cluster, with all associated mutation names
  cluster_labels <- variants %>%
    group_by(clone_id, original_cluster) %>%
    summarise(
      mutations = paste(unique(mutation_label), collapse = "; "),
      .groups = "drop"
    )
  
  list(variants = variants, cluster_labels = cluster_labels)
}


# Label the Jellyfisher phylogeny with mutations

label_phylogeny <- function(tables, cluster_labels) {
  
  phy_ids <- as.character(tables$phylogeny$subclone)
  
  # Jellyfisher uses subclone IDs to connect its tables.
  # Rename IDs consistently in phylogeny and compositions.
  label_map <- setNames(phy_ids, phy_ids)
  label_map["0"] <- "Root (0)"
  
  for (i in seq_len(nrow(cluster_labels))) {
    
    id <- as.character(cluster_labels$clone_id[i])
    
    if (id %in% phy_ids) {
      label_map[id] <- paste0(
        "Clone ", id,
        " [SciClone cluster ", cluster_labels$original_cluster[i],
        "]: ",
        cluster_labels$mutations[i]
      )
    }
  }
  
  relabel <- function(x) {
    x <- as.character(x)
    result <- unname(label_map[x])
    result[is.na(result)] <- x[is.na(result)]
    result
  }
  
  tables$phylogeny$subclone <- relabel(
    tables$phylogeny$subclone
  )
  
  tables$phylogeny$parent <- relabel(
    tables$phylogeny$parent
  )
  
  tables$compositions$subclone <- relabel(
    tables$compositions$subclone
  )
  
  tables
}


# Process every saved ClonEvol result
rds_files <- list.files(
  "./clonevol_results/objects/",
  pattern = "_clonevol\\.rds$",
  full.names = TRUE
)

if (length(rds_files) == 0) {
  stop("No ClonEvol RDS files found")
}

process_patient <- function(rds_file) {
  
  patient_id <- sub(
    "_clonevol\\.rds$", "",
    basename(rds_file)
  )
  
  message("\nProcessing ", patient_id)
  
  tryCatch({
    
    y <- readRDS(rds_file)
    
    # Only plot patients with an available consensus tree
    n_trees <- length(y$matched$merged.trees)
    
    if (n_trees == 0) {
      message("Skipping: no matched consensus tree.")
      
      return(data.frame(
        Patient_ID = patient_id,
        Status = "Skipped: no consensus tree",
        Error = NA_character_
      ))
    }
    
    # Convert the tree and extract Jellyfisher tables
    y_branch <- clonevol::convert.consensus.tree.clone.to.branch(y)
    
    tables <- jellyfisher::extract_tables_from_clonevol(
      y_branch,
      model = 1
    )
    
    # Match actual mutations to the inferred clone IDs
    mutation_data <- prepare_mutation_labels(patient_id)
    
    tables <- label_phylogeny(
      tables,
      mutation_data$cluster_labels
    )
    
    # Save the mutation-to-cluster lookup for verification
    write_tsv(
      mutation_data$cluster_labels,
      file.path("./jellyfisher_results/tables/", paste0(patient_id, "_cluster_mutations.tsv"))
    )
    
    # Create and save the interactive plot
    widget <- jellyfisher::jellyfisher(
      tables,
      height = "700px",
      width = "100%"
    )
    
    # Temporary HTML file needed to render the widget
    temp_html <- file.path("./jellyfisher_results/plots/",
                           paste0(patient_id,"_temp.html"))
    
    # Final PNG image
    png_file <- file.path(
      "./jellyfisher_results/plots/",
      paste0(patient_id, "_jellyfisher.png")
    )
    
    # Save the widget temporarily as HTML
    htmlwidgets::saveWidget(
      widget,
      file = temp_html,
      selfcontained = TRUE
    )
    
    # Capture the rendered widget as a PNG
    webshot2::webshot(
      url = temp_html,
      file = png_file,
      vwidth = 1600,
      vheight = 1000,
      delay = 3
    )
    
    # Remove the temporary HTML file after successful capture
    file.remove(temp_html)
    
    if (show_in_viewer) {
      print(widget)
    }
    
    message("Saved Jellyfisher plot for ", patient_id)
    
    data.frame(
      Patient_ID = patient_id,
      Status = "Plot saved",
      Error = NA_character_
    )
    
  }, error = function(e) {
    
    message("Failed for ", patient_id, ": ", conditionMessage(e))
    
    data.frame(
      Patient_ID = patient_id,
      Status = "Failed",
      Error = conditionMessage(e)
    )
  })
}

# Run all saved patient objects
plot_summary <- bind_rows(lapply(rds_files, process_patient))

write_tsv(
  plot_summary,
  file.path("./jellyfisher_results/", "jellyfisher_run_summary.tsv")
)

# Final results
print(plot_summary)
print(count(plot_summary, Status))

