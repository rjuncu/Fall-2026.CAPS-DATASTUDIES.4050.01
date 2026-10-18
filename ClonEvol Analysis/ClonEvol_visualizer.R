# install.packages("jellyfisher")
# install.packages("devtools")
# devtools::install_github("hdng/clonevol")
# install.packages("gridBase")
# install.packages("gridExtra")
# install.packages("ggplot2")
# install.packages("igraph")
# install.packages("packcircles")
# devtools::install_github("hdng/trees")
# install.packages("readr")
# install.packages("dplyr")
# install.packages("htmlwidgets")
# install.packages("webshot2")
# install.packages("magick")

library(readr)
library(dplyr)
library(clonevol)
library(jellyfisher)
library(htmlwidgets)
library(webshot2)
library(magick)

# Create output directories
dir.create(
  "jellyfisher_results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "jellyfisher_results/plots",
  recursive = TRUE,
  showWarnings = FALSE
)

# Set TRUE to open each interactive plot in the RStudio Viewer
show_in_viewer <- FALSE

# Input files
sciclone_total <- read_tsv(
  "C:/Users/rjunc/Documents/4050_project/Fall-2026.CAPS-DATASTUDIES.4050.01/SciClone Analysis/sciclone_results/tables/sciclone_total.tsv",
  show_col_types = FALSE
)

annotations <- read_tsv(
  "C:/Users/rjunc/Documents/4050_project/Fall-2026.CAPS-DATASTUDIES.4050.01/SciClone Analysis/annotated_mutations.tsv",
  show_col_types = FALSE
)

# Check required columns
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

# Normalize chromosome and position values for joining
normalize_chr <- function(x) {
  toupper(
    sub(
      "^chr", "",
      trimws(as.character(x)),
      ignore.case = TRUE
    )
  )
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
  c(
    "HGVSp_Short", "HGVSp", "Protein_Change",
    "ProteinChange", "AAChange", "Amino_acids"
  ),
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

message(
  "Protein-change column: ",
  ifelse(is.na(protein_col), "not found", protein_col)
)

# Create annotation lookup
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
      z <- unique(
        mutation_label[
          !is.na(mutation_label) & mutation_label != ""
        ]
      )
      
      if (length(z)) {
        paste(z, collapse = " / ")
      } else {
        NA_character_
      }
    },
    .groups = "drop"
  )

# Prepare mutation annotations for one patient
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
  
  # Create consecutive SciClone cluster labels for reference.
  # These are NOT assumed to be ClonEvol branch IDs.
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
  
  # One row per SciClone cluster
  cluster_labels <- variants %>%
    group_by(clone_id, original_cluster) %>%
    summarise(
      mutations = paste(unique(mutation_label), collapse = "; "),
      .groups = "drop"
    )
  
  list(
    variants = variants,
    cluster_labels = cluster_labels
  )
}

# Locate saved ClonEvol results
rds_files <- list.files(
  "./clonevol_results/objects/",
  pattern = "_clonevol\\.rds$",
  full.names = TRUE
)

if (length(rds_files) == 0) {
  stop("No ClonEvol RDS files found.")
}

# Process one patient
process_patient <- function(rds_file) {
  
  patient_id <- sub(
    "_clonevol\\.rds$", "",
    basename(rds_file)
  )
  
  message("\nProcessing ", patient_id)
  
  # Track temporary files so they can be cleaned up on failure too
  temp_html <- file.path(
    "./jellyfisher_results/plots/",
    paste0(patient_id, "_temp.html")
  )
  
  tree_png <- file.path(
    "./jellyfisher_results/plots/",
    paste0(patient_id, "_tree_temp.png")
  )
  
  tryCatch({
    
    y <- readRDS(rds_file)
    
    # Check for a matched consensus tree
    n_trees <- length(y$matched$merged.trees)
    
    if (n_trees == 0) {
      message("Skipping: no matched consensus tree.")
      
      return(data.frame(
        Patient_ID = patient_id,
        Status = "Skipped: no consensus tree",
        Error = NA_character_
      ))
    }
    
    # Convert tree and extract Jellyfisher tables
    y_branch <- clonevol::convert.consensus.tree.clone.to.branch(y)
    
    tables <- jellyfisher::extract_tables_from_clonevol(
      y_branch,
      model = 1
    )
    
    # Prepare mutation annotations
    mutation_data <- prepare_mutation_labels(patient_id)
    
    # Save mutation-to-SciClone-cluster lookup
    write_tsv(
      mutation_data$cluster_labels,
      file.path(
        "./jellyfisher_results/tables/",
        paste0(patient_id, "_cluster_mutations.tsv")
      )
    )
    
    # Create Jellyfisher widget using the original clone IDs
    widget <- jellyfisher::jellyfisher(
      tables,
      height = "700px",
      width = "100%"
    )
    
    # Capture the widget as a temporary PNG
    htmlwidgets::saveWidget(
      widget,
      file = temp_html,
      selfcontained = TRUE
    )
    
    webshot2::webshot(
      url = temp_html,
      file = tree_png,
      vwidth = 1400,
      vheight = 1000,
      delay = 3
    )
    
    # Build the mutation legend
    legend_data <- mutation_data$cluster_labels %>%
      arrange(clone_id)
    
    legend_lines <- unlist(
      lapply(seq_len(nrow(legend_data)), function(i) {
        
        mutations <- strwrap(
          legend_data$mutations[i],
          width = 42
        )
        
        c(
          paste0(
            "SciClone cluster ",
            legend_data$original_cluster[i]
          ),
          paste0("  - ", mutations),
          ""
        )
      }),
      use.names = FALSE
    )
    
    legend_text <- paste(
      c(
        paste("Patient:", patient_id),
        "MUTATION LEGEND",
        "==============================",
        "",
        legend_lines
      ),
      collapse = "\n"
    )
    
    # Estimate a generous height for the legend
    legend_height <- max(
      1200,
      45 * (length(strsplit(legend_text, "\n")[[1]]) + 2)
    )
    
    # Create legend panel
    legend_img <- magick::image_blank(
      width = 650,
      height = legend_height,
      color = "white"
    )
    
    legend_img <- magick::image_annotate(
      legend_img,
      text = legend_text,
      gravity = "northwest",
      location = "+25+25",
      size = 18,
      font = "Arial",
      color = "black"
    )
    
    # Read tree screenshot
    tree_img <- magick::image_read(tree_png)
    
    # Extend tree canvas to match legend height
    tree_width <- magick::image_info(tree_img)$width
    
    tree_img <- magick::image_extent(
      tree_img,
      geometry = paste0(tree_width, "x", legend_height),
      gravity = "north",
      color = "white"
    )
    
    # Combine tree and legend side by side
    combined_img <- magick::image_append(
      c(tree_img, legend_img),
      stack = FALSE
    )
    
    # Save final PNG
    final_png <- file.path(
      "./jellyfisher_results/plots/",
      paste0(patient_id, "_jellyfisher.png")
    )
    
    magick::image_write(
      combined_img,
      path = final_png,
      format = "png"
    )
    
    # Clean up intermediate files
    file.remove(temp_html, tree_png)
    
    message("Saved tree with mutation legend: ", final_png)
    
    if (show_in_viewer) {
      print(widget)
    }
    
    data.frame(
      Patient_ID = patient_id,
      Status = "Plot saved",
      Error = NA_character_
    )
    
  }, error = function(e) {
    
    # Remove temporary files if processing fails
    if (file.exists(temp_html)) {
      file.remove(temp_html)
    }
    
    if (file.exists(tree_png)) {
      file.remove(tree_png)
    }
    
    message(
      "Failed for ", patient_id, ": ",
      conditionMessage(e)
    )
    
    data.frame(
      Patient_ID = patient_id,
      Status = "Failed",
      Error = conditionMessage(e)
    )
  })
}

# Run all patient objects
plot_summary <- bind_rows(
  lapply(rds_files, process_patient)
)

# Save batch summary
write_tsv(
  plot_summary,
  "./jellyfisher_results/jellyfisher_run_summary.tsv"
)

# Display results
print(plot_summary)
print(count(plot_summary, Status))
