# clonevol_cluster_analysis.R  (v4 - console only; nothing is written to disk)
#
# Cluster proportions + hugo per cluster + bootstrap tree position + kinetics
# bands + MSigDB functional annotation + cross-patient crosstabs + the
# mutation-sharing premise check.
#
# RUN FROM THE REPOSITORY ROOT (the folder containing "SciClone Analysis" and
# "ClonEvol Analysis"); the guard below stops a wrong getwd() loudly. No file is
# created, modified or deleted: every table is printed AND left in the global
# environment under the name listed in the footer, so you export what you want.
#
# DESIGN: kinetics (A), function (B) and tree position (C) are classified
# separately and only then crossed - asserting the phase order inside the rubric
# would only recover the assumption instead of testing it.

ensure_pkg <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressMessages(library(p, character.only = TRUE))
}
invisible(lapply(c("readr", "dplyr", "tidyr", "msigdbr"), ensure_pkg))

sect <- function(...) cat("\n== ", ..., " ==\n")

sciclone_file  <- "./SciClone Analysis/sciclone_results/tables/sciclone_total.tsv"
annotated_file <- "./SciClone Analysis/annotated_mutations.tsv"
objects_dir    <- "./ClonEvol Analysis/clonevol_results/objects"
for (p in c(sciclone_file, annotated_file, objects_dir)) {
  if (!file.exists(p))
    stop("'", p, "' not found from ", getwd(), ". Run from the repository root.")
}


DEPTH_MIN <- 100   # a variant needs >= 100 reads before "absent" means anything

maf <- read_tsv("./SciClone Analysis/annotated_mutations.tsv", show_col_types = FALSE) %>%
  filter((t_ref_count + t_alt_count) >= DEPTH_MIN) %>%
  mutate(patient     = sub("-[PM]$", "", Tumor_Sample_Barcode),   # EV-001-P -> EV-001
         Mutation_ID = paste(Chromosome, Start_Position, sep = ":"))

sites <- maf %>%
  group_by(patient, Mutation_ID) %>%
  summarise(in_primary    = any(`Sample Type` == "Primary"),
            in_metastasis = any(`Sample Type` == "Metastasis"),
            hugo          = paste(sort(unique(Hugo_Symbol)), collapse = "; "),
            .groups = "drop") %>%
  mutate(sharing = case_when(in_primary & in_metastasis ~ "shared",
                             in_primary                ~ "primary_only",
                             TRUE                      ~ "metastasis_only"))

met_only <- sites %>% filter(sharing == "metastasis_only")

# Same sites, but exploded back to one row per observed mutation so you can see
# depth and VAF in the metastasis sample.
met_only_rows <- maf %>%
  inner_join(met_only %>% select(patient, Mutation_ID), by = c("patient", "Mutation_ID")) %>%
  transmute(patient, Mutation_ID, hugo = Hugo_Symbol, chr = Chromosome, st = Start_Position,
            sample = Tumor_Sample_Barcode, alt = t_alt_count,
            depth = t_ref_count + t_alt_count, vaf = t_alt_count / (t_ref_count + t_alt_count))

print(table(sites$sharing))
print(as.data.frame(met_only))

{
  gc  <- head(sort(table(met_only$hugo), decreasing = TRUE), 5)
  pos <- barplot(as.numeric(gc), horiz = TRUE, las = 1, cex.names = 0.9,
                 names.arg = names(gc), col = "#6b8fbd", border = "white",
                 xlim = c(0, max(gc) + 0.8),
                 main = "Top 5 genes by metastasis-only mutation sites",
                 xlab = "mutation sites absent from that patient's primary tumor")
  text(as.numeric(gc) + 0.12, pos, labels = as.character(gc), cex = 0.9)
}


# Coerce a hit/flag column: missing -> 0, NA -> 0, single-cell list column ->
# flattened, multi-value cell -> stop (that means a join duplicated the row, so
# every count built on it would be inflated). Used by Parts 4 and 5.
flag_cols <- function(df, cols) {
  for (cc in cols) {
    v <- if (cc %in% names(df)) df[[cc]] else 0L
    if (is.list(v)) {
      if (any(lengths(v) != 1L))
        stop("'", cc, "': one row carries more than one value - a join ",
             "duplicated rows and the counts would be wrong.")
      v <- unlist(v)
    }
    df[[cc]] <- as.integer(ifelse(is.na(v), 0L, v))
  }
  df
}

# ----------------------------------------------------------------- rubric ----
# Axis A - kinetics: ratio = mean cluster VAF (met) / mean cluster VAF (primary).
# 0.67 = 1/1.5, so both directions weigh the same on the log scale. Fixed in
# advance, never tuned to the results.
RATIO_PRIMARY_ENRICHED <- 0.67
RATIO_MET_ENRICHED     <- 1.5

# Axis B - function: MSigDB gene sets matched by NAME pattern. What actually
# matched lands in `matched_sets`, so the rubric stays auditable.
CATEGORY_PATTERNS <- list(
  ADHESION_ECM   = c("INTEGRIN", "CELL_SURFACE_INTERACTION", "CELL_ADHESION",
                     "FOCAL_ADHESION", "APICAL_JUNCTION", "APICAL_SURFACE",
                     "EXTRACELLULAR_MATRIX", "ECM"),
  PROLIFERATION  = c("CELL_CYCLE", "MITOTIC", "E2F_TARGETS", "MYC_TARGETS",
                     "G2M_CHECKPOINT", "MTORC1_SIGNALING", "DNA_REPLICATION"),
  WNT_INITIATION = c("WNT_B_CATENIN", "WNT_SIGNALING", "WNT")
)

# Axis C - tree position: founder-like if founding in >= 50% of that patient's
# bootstrap models, subclone-like if subclonal in >= 50%, otherwise unstable.
FOUNDER_FRAC_CUTOFF <- 0.5

# Tree labs are matched to clusters by VAF, not by index (see Part 4). A match
# has to fall within this many percentage points of the cluster mean primary VAF.
TREE_VAF_TOL <- 1.0

# Sanity anchor: EV-001 must show 50.90406, 20.19231, 20.10856, 15.72519.
anchor <- file.path(objects_dir, "EV-001_clonevol.rds")
if (file.exists(anchor)) {
  ev_001 <- readRDS(anchor)
  cat("EV-001 bootstrap models:", length(ev_001$models$Primary),
      "| matched models:", ev_001$num.matched.models, "\n")
}

# =========================================== 1. eligible variants + hugo =====
# Eligible = the rows clonevol actually clustered, renumbered the way clonevol
# renumbered them. Labels are PATIENT-LOCAL: never compare a lab across patients.
eligible <- read_tsv(sciclone_file) %>%
  filter(adequateDepth == 1, !is.na(cluster), cluster > 0,
         !is.na(Primary_VAF), !is.na(Metastasis_VAF)) %>%
  distinct(Patient_ID, chr, st, .keep_all = TRUE) %>%
  mutate(lab = as.character(match(cluster, sort(unique(cluster))) - 1L))

# sciclone_total.tsv carries the bare Patient_ID (EV-001) while
# annotated_mutations.tsv uses EV-001-P / EV-001-M barcodes: join the primary
# barcode, then fall back to the metastasis one.
genes <- read_tsv(annotated_file) %>%
  transmute(sample = Tumor_Sample_Barcode, chr = as.character(Chromosome),
            st = as.numeric(Start_Position), hugo = Hugo_Symbol)

variants <- eligible %>%
  mutate(barcode_p = paste0(Patient_ID, "-P"), barcode_m = paste0(Patient_ID, "-M")) %>%
  left_join(genes, by = c("barcode_p" = "sample", "chr" = "chr", "st" = "st")) %>%
  left_join(genes, by = c("barcode_m" = "sample", "chr" = "chr", "st" = "st"),
            suffix = c("", "_m")) %>%
  mutate(hugo = coalesce(hugo, hugo_m)) %>%
  select(-hugo_m, -barcode_p, -barcode_m)

cat("Eligible variants:", nrow(variants), "| no hugo match:", sum(is.na(variants$hugo)), "\n")
# $variants inside each .rds is NOT used for cluster->gene mapping: its cluster
# column follows the bootstrap resample, not the original labels.

# ================================= 2. clusters, proportions, kinetics ========
# 2a. One row per patient-cluster: variant count, gene count, gene list, mean VAFs
clusters <- variants %>%
  group_by(Patient_ID, lab) %>%
  summarise(n_variants = n(),
            n_genes    = dplyr::n_distinct(hugo),
            gene_list  = paste(sort(unique(na.omit(hugo))), collapse = "; "),
            mean_vaf_primary    = mean(Primary_VAF),
            mean_vaf_metastasis = mean(Metastasis_VAF),
            .groups = "drop")
# 2b. Each cluster's share of its own tumour's total VAF
clusters <- clusters %>%
  group_by(Patient_ID) %>%
  mutate(prop_primary    = mean_vaf_primary    / sum(mean_vaf_primary),
         prop_metastasis = mean_vaf_metastasis / sum(mean_vaf_metastasis)) %>%
  ungroup()
# 2c. Kinetics = metastasis mean VAF divided by primary mean VAF.
# A primary mean of 0 becomes NA, so no divide-by-zero warnings and no Inf.
primary_mean <- ifelse(clusters$mean_vaf_primary > 0, clusters$mean_vaf_primary, NA_real_)
clusters$ratio_met_over_primary <- clusters$mean_vaf_metastasis / primary_mean
# 2d. Three bands, fixed before looking at the data
clusters$kinetic_class <-
  ifelse(is.na(clusters$ratio_met_over_primary),                      "unclassifiable_kinetics",
         ifelse(clusters$ratio_met_over_primary < RATIO_PRIMARY_ENRICHED,   "primary_enriched",
                ifelse(clusters$ratio_met_over_primary > RATIO_MET_ENRICHED,       "metastasis_enriched",
                       "shared_similar")))
clusters <- clusters %>% arrange(Patient_ID, as.integer(lab))
# 2e. The two tables you export from
cluster_proportions <- clusters %>%
  select(Patient_ID, lab, n_variants, mean_vaf_primary, mean_vaf_metastasis,
         prop_primary, prop_metastasis, ratio_met_over_primary, kinetic_class)
cluster_genes <- clusters %>%
  select(Patient_ID, lab, n_variants, n_genes, gene_list)
# 2f. Band table and band chart, same bins, top bin left open so nothing is dropped
sect("kinetic band sensitivity (0.67 and 1.50 = shared_similar)")
bands       <- cut(clusters$ratio_met_over_primary,
                   breaks = c(0, 0.5, 0.67, 1, 1.5, 2, Inf), right = FALSE,
                   labels = c("<0.50", "0.50-0.67", "0.67-1.00", "1.00-1.50", "1.50-2.00", ">=2"))
band_counts <- table(bands)
print(band_counts)
pos <- barplot(as.numeric(band_counts), col = "#6b8fbd", border = "white", space = 0.15,
               names.arg = names(band_counts), cex.names = 0.85,
               ylim = c(0, max(band_counts) * 1.18),
               xlab = "cluster VAF ratio band (metastasis mean / primary mean)",
               ylab = "number of clusters", main = "Kinetic band distribution")
text(pos, as.numeric(band_counts) + 2, labels = as.character(band_counts), cex = 0.95)


# ====================================== 3. MSigDB functional annotation ======
fetch_sets <- function(collection) {
  d <- tryCatch(msigdbr::msigdbr(collection = collection),
                error = function(e)
                  tryCatch(msigdbr::msigdbr(category = collection),
                           error = function(e2) NULL))
  if (is.null(d))
    stop("Could not fetch MSigDB collection '", collection,
         "'. Check the connection, or install.packages(\"msigdbr\").")
  transmute(d, gs_name = toupper(gs_name), gene = .data[[
    if ("gene_symbol" %in% names(d)) "gene_symbol" else "human_gene_symbol"]])
}
all_sets <- bind_rows(fetch_sets("H"), fetch_sets("C2")) %>% distinct(gs_name, gene)

matched_sets <- bind_rows(lapply(names(CATEGORY_PATTERNS), function(cat)
  data.frame(category = cat, stringsAsFactors = FALSE,
             gs_name = unique(unlist(lapply(CATEGORY_PATTERNS[[cat]], grep,
                                            unique(all_sets$gs_name), value = TRUE))))))
sect("gene sets matched per category"); print(table(matched_sets$category))

present_genes <- sort(unique(variants$hugo[!is.na(variants$hugo)]))
if (length(present_genes) == 0)
  stop("No hugo symbols matched any variant - check the Hugo_Symbol column of annotated_mutations.tsv.")

gene_category <- matched_sets %>% inner_join(all_sets, by = "gs_name") %>%
  filter(gene %in% present_genes) %>% distinct(category, gene) %>% mutate(hit = TRUE)

if (nrow(gene_category) == 0)
  cat("WARNING: no observed gene matched any gene set; every cluster will be",
      "UNANNOTATED and the tests will be skipped.\n")

# Mirror audit: every observed gene, the rejected ones included, with the reason.
gene_audit <- expand.grid(gene = present_genes, category = names(CATEGORY_PATTERNS),
                          stringsAsFactors = FALSE) %>%
  left_join(gene_category, by = c("gene", "category")) %>%
  group_by(gene) %>%
  summarise(categories = paste(category[!is.na(hit)], collapse = "; "), .groups = "drop") %>%
  mutate(status = ifelse(categories == "", "excluded: not in any matched gene set",
                         ifelse(grepl("; ", categories),
                                "ambiguous: belongs to >1 category", "annotated")))
cat("Genes observed:", nrow(gene_audit), "| annotated:", sum(gene_audit$status == "annotated"),
    "| ambiguous:", sum(grepl("ambiguous", gene_audit$status)),
    "| excluded:", sum(grepl("excluded", gene_audit$status)), "\n")

{
  genes_per_cat <- table(gene_category$category)
  audit_status  <- table(factor(gene_audit$status,
                                levels = c("annotated", "ambiguous: belongs to >1 category",
                                           "excluded: not in any matched gene set")))
  
  par(mfcol = c(1, 2), mar = c(5, 9, 3.5, 1))
  
  p1 <- barplot(as.numeric(genes_per_cat), horiz = TRUE, las = 1, col = "#6b8fbd",
                border = "white", xlim = c(0, max(genes_per_cat) * 1.25),
                names.arg = names(genes_per_cat), cex.names = 0.8,
                main = "Observed genes per functional category", xlab = "genes")
  text(as.numeric(genes_per_cat) + 1, p1, labels = as.character(genes_per_cat), cex = 0.95)
  
  p2 <- barplot(as.numeric(audit_status), horiz = TRUE, las = 1,
                col = c("#6b8fbd", "#d9a066", "#b25b4c"), border = "white",
                xlim = c(0, max(audit_status) * 1.25), names.arg = c("annotated", "ambiguous", "excluded"),
                cex.names = 0.8, main = "Audit of all observed genes", xlab = "genes")
  text(as.numeric(audit_status) + 1, p2, labels = as.character(audit_status), cex = 0.95)
  
  mtext("ambiguous genes are counted in every category they match", side = 1,
        outer = TRUE, line = 2.6, cex = 0.8)
}


# Paste this AFTER the MSigDB section of your script. Needs three objects from
# there: variants (with hugo + lab), gene_category (category, gene), gene_audit.
# Prints the table and leaves it in `top5`. Writes no files.

missing <- setdiff(c("variants", "gene_category", "gene_audit"), ls())
if (length(missing))
  stop("Paste this after the sections that create: ", paste(missing, collapse = ", "))

ambiguous_genes <- gene_audit$gene[grepl("ambiguous", gene_audit$status)]

# One row per (category, gene). An ambiguous gene appears once under each
# category it belongs to, which is what the cluster hit counts already did.
cat_gene <- variants %>%
  filter(!is.na(hugo)) %>%
  inner_join(gene_category %>% select(category, gene),
             by = c("hugo" = "gene"), relationship = "many-to-many") %>%
  group_by(category, hugo) %>%
  summarise(n_variants = n(),                                          # mutation rows
            n_sites    = dplyr::n_distinct(paste(chr, st, sep = ":")), # distinct positions
            n_patients = dplyr::n_distinct(Patient_ID),
            n_clusters = dplyr::n_distinct(paste(Patient_ID, lab, sep = "/")),
            .groups = "drop") %>%
  mutate(ambiguous = hugo %in% ambiguous_genes)

# Category sizes measured on the FULL table, before slicing down to 5.
cat_totals <- cat_gene %>%
  group_by(category) %>%
  summarise(genes_in_category = dplyr::n_distinct(hugo),
            genes_single_category = sum(!(hugo %in% ambiguous_genes)), .groups = "drop")

# Rank by distinct genomic positions; ties broken by patients, then alphabetical.
top5 <- cat_gene %>%
  group_by(category) %>%
  slice_max(order_by = n_sites, n = 5, with_ties = FALSE) %>%
  arrange(category, desc(n_sites), desc(n_patients), hugo) %>%
  group_by(category) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  left_join(cat_totals, by = "category") %>%
  select(category, rank, hugo, n_variants, n_sites, n_patients, n_clusters,
         ambiguous, genes_in_category, genes_single_category)

cat("\n==  TOP 5 GENES PER CATEGORY (ranked by distinct genomic positions)  ==\n")
print(as.data.frame(top5))


# ======================================= 4. bootstrap tree position per lab ==
model_counts <- list(); tree_meta <- list(); trees <- list()
for (f in list.files(objects_dir, pattern = "\\.rds$", full.names = TRUE)) {
  pid <- sub("_clonevol\\.rds$", "", basename(f))
  y <- tryCatch(readRDS(f), error = function(e) NULL)        # 43-byte NULL objects
  m <- if (is.null(y)) NULL else y$models[["Primary"]]
  model_counts[[pid]] <- data.frame(Patient_ID = pid, n_bootstrap_models = length(m),
                                    num_matched_models = if (is.null(y)) NA_integer_ else y$num.matched.models)
  if (length(m) == 0) next
  
  tree_meta[[pid]] <- bind_rows(lapply(m, function(tr) tr[, c("lab", "is.founder", "is.subclone", "vaf")])) %>%
    mutate(lab = as.character(lab), Patient_ID = pid) %>%
    group_by(Patient_ID, lab) %>%
    summarise(n_models_with_lab = n(), founder_frac = mean(is.founder),
              subclone_frac = mean(is.subclone), tree_vaf = mean(vaf), .groups = "drop") %>%
    mutate(tree_role = case_when(founder_frac >= FOUNDER_FRAC_CUTOFF  ~ "founder-like",
                                 subclone_frac >= FOUNDER_FRAC_CUTOFF ~ "subclone-like",
                                 TRUE                                  ~ "unstable_topology"))
  trees[[pid]] <- cbind(m[[1]][, c("lab", "parent", "ancestors", "vaf", "is.founder", "is.subclone")],
                        Patient_ID = pid) %>% mutate(lab = as.character(lab))
}
model_counts <- bind_rows(model_counts); tree_meta <- bind_rows(tree_meta); trees <- bind_rows(trees)
cat("Patients with NULL/empty objects:", sum(model_counts$n_bootstrap_models == 0),
    "| clusters with tree metadata:", nrow(tree_meta), "\n")

# clonevol clustered MORE clusters than the eligibility filter keeps. EV-001's
# bootstrap trees have five labs (50.90 / 25.59 / 20.19 / 20.11 / 15.73) while
# the filtered table keeps four (25.6 / 20.2 / 20.1 / 15.7), so renumbering
# shifts by one and a tree lab must be matched to a cluster by VAF - never by
# its index. Joining on lab would attach the wrong founder/subclone history.
tree_meta <- tree_meta %>% mutate(tree_lab = lab, tree_vaf_pct = tree_vaf * 100)

clusters <- bind_rows(lapply(split(clusters, clusters$Patient_ID), function(cl) {
  tm <- tree_meta %>% filter(Patient_ID == cl$Patient_ID[1])
  if (nrow(tm) == 0)
    return(mutate(cl, tree_lab = NA_character_, n_models_with_lab = NA_integer_,
                  founder_frac = NA_real_, subclone_frac = NA_real_,
                  tree_role = "no_consensus_models"))
  k <- vapply(cl$mean_vaf_primary, function(v) {
    d <- abs(tm$tree_vaf_pct - v); i <- which.min(d)
    if (d[i] > TREE_VAF_TOL) NA_integer_ else i
  }, integer(1))
  cl$tree_lab          <- tm$tree_lab[k]
  cl$n_models_with_lab <- tm$n_models_with_lab[k]
  cl$founder_frac      <- tm$founder_frac[k]
  cl$subclone_frac     <- tm$subclone_frac[k]
  cl$tree_role         <- ifelse(is.na(k), "tree_vaf_unmatched", tm$tree_role[k])
  cl
})) %>% arrange(Patient_ID, as.integer(lab))

sect("tree lab alignment")
print(table(no_tree = clusters$tree_role == "no_consensus_models",
            unmatched = clusters$tree_role == "tree_vaf_unmatched",
            matched = clusters$tree_role %in% c("founder-like", "subclone-like", "unstable_topology")))
sect("EV-001 cluster -> tree lab mapping")
print(clusters %>% filter(Patient_ID == "EV-001") %>%
        select(lab, mean_vaf_primary, tree_lab, founder_frac, subclone_frac, tree_role))

# ======================================== 5. cross the three axes ============
hits <- variants %>%
  distinct(Patient_ID, lab, hugo) %>% filter(!is.na(hugo)) %>%
  inner_join(gene_category, by = c("hugo" = "gene"), relationship = "many-to-many") %>%
  group_by(Patient_ID, lab, category) %>%
  summarise(n = dplyr::n_distinct(hugo), .groups = "drop")
hits <- if (nrow(hits) > 0) {
  pivot_wider(hits, names_from = category, values_from = n, values_fill = 0, names_prefix = "n_")
} else {
  data.frame(Patient_ID = character(0), lab = character(0))
}

cat_cols <- paste0("n_", names(CATEGORY_PATTERNS))
clusters <- flag_cols(clusters %>% left_join(hits, by = c("Patient_ID", "lab")), cat_cols)

hit_mat <- as.matrix(clusters[, cat_cols])
clusters$total_hits <- rowSums(hit_mat)
clusters$n_tied     <- rowSums(hit_mat == do.call(pmax, as.data.frame(hit_mat)))
# Empty and tied rows are decided BEFORE any argmax: which.max would stamp them
# with whichever category happens to be listed first.
clusters$functional_class <- ifelse(clusters$total_hits == 0, "UNANNOTATED",
                                    ifelse(clusters$n_tied > 1, "MIXED_AMBIGUOUS", names(CATEGORY_PATTERNS)[max.col(hit_mat)]))

profile <- clusters %>%
  select(Patient_ID, lab, n_variants, n_genes, gene_list, mean_vaf_primary,
         mean_vaf_metastasis, prop_primary, prop_metastasis, ratio_met_over_primary,
         kinetic_class, tree_lab, n_models_with_lab, founder_frac, subclone_frac, tree_role,
         all_of(cat_cols), total_hits, n_tied, functional_class)

# Mirror exclusions: every cluster that was not placed, with its reason.
excluded <- profile %>%
  filter(functional_class %in% c("UNANNOTATED", "MIXED_AMBIGUOUS")) %>%
  mutate(reason = ifelse(functional_class == "UNANNOTATED",
                         "no cluster gene matched any matched gene set (likely passengers)",
                         "categories tied on hit count - deliberately unassigned")) %>%
  select(Patient_ID, lab, gene_list, n_variants, n_genes, all_of(cat_cols), reason)

classified <- profile %>% filter(!functional_class %in% c("UNANNOTATED", "MIXED_AMBIGUOUS"))
sect("classification counts")
cat("clusters:", nrow(profile), "| classified:", nrow(classified), "| excluded:", nrow(excluded), "\n")

sect("functional class x kinetics")
print(table(functional_class = classified$functional_class, kinetic_class = classified$kinetic_class))
sect("functional class x tree position")
print(table(functional_class = classified$functional_class, tree_role = classified$tree_role))

two_phase <- classified %>% filter(functional_class %in% c("ADHESION_ECM", "PROLIFERATION"))
if (nrow(two_phase) >= 4) {
  ct <- table(two_phase$functional_class, two_phase$kinetic_class == "metastasis_enriched")
  sect("Fisher: adhesion vs proliferation x metastasis-enriched or not")
  print(ct); print(fisher.test(ct))
}
if (length(unique(classified$functional_class)) > 1) {
  sect("Kruskal: met/primary VAF ratio by functional class")
  print(kruskal.test(ratio_met_over_primary ~ functional_class, data = classified))
  founder_data <- classified %>% filter(!is.na(founder_frac))
  if (length(unique(founder_data$functional_class)) > 1) {
    sect("Kruskal: founder_frac by functional class")
    print(kruskal.test(founder_frac ~ functional_class, data = founder_data))
  }
}

# ================================ 6. premise check on the full mutation table
# "Metastasis mutations differ from primary" is NOT testable from the cluster
# objects: adequateDepth == 1 required coverage in BOTH samples, so metastasis-only
# variants were deleted before clustering. This tests the premise directly.
mut_all <- read_tsv(annotated_file) %>%
  filter((t_ref_count + t_alt_count) >= 100) %>%
  mutate(patient = sub("-[PM]$", "", Tumor_Sample_Barcode),
         Mutation_ID = paste(Chromosome, Start_Position, sep = ":"))

# Category flags per site, then one row per mutation SITE: a site counts once per
# category it touches, and genes/annotations can never multiply the unit counted.
site_flags <- {
  sc <- distinct(mut_all, patient, Mutation_ID, hugo = Hugo_Symbol) %>%
    inner_join(gene_category, by = c("hugo" = "gene"), relationship = "many-to-many") %>%
    distinct(patient, Mutation_ID, category)
  if (nrow(sc) > 0) {
    mutate(sc, touch = 1) %>%
      pivot_wider(names_from = category, values_from = touch, values_fill = 0,
                  names_prefix = "cat_")
  } else {
    data.frame(patient = character(0), Mutation_ID = character(0))
  }
}

mut_sites <- mut_all %>%
  group_by(patient, Mutation_ID) %>%
  summarise(in_primary = any(`Sample Type` == "Primary"),
            in_metastasis = any(`Sample Type` == "Metastasis"), .groups = "drop") %>%
  mutate(sharing = case_when(in_primary & in_metastasis ~ "shared",
                             in_primary ~ "primary_only",
                             TRUE ~ "metastasis_only")) %>%
  left_join(site_flags, by = c("patient", "Mutation_ID"))

cat_cols2 <- paste0("cat_", names(CATEGORY_PATTERNS))
mut_sites <- flag_cols(mut_sites, cat_cols2)

sect("sharing of mutation sites"); print(table(mut_sites$sharing))
cat("total sites:", nrow(mut_sites), "| annotated sites (>=1 category):",
    sum(rowSums(as.matrix(mut_sites[, cat_cols2])) > 0), "\n")

ann <- mut_sites %>%
  pivot_longer(all_of(cat_cols2), names_to = "category", values_to = "touch") %>%
  filter(touch == 1) %>% mutate(category = sub("^cat_", "", category))
sect("category x sharing (a site counts once per category it touches)")
print(table(category = ann$category, sharing = ann$sharing))

mp <- ann %>% filter(category %in% c("ADHESION_ECM", "PROLIFERATION"),
                     sharing %in% c("shared", "metastasis_only"))
if (nrow(mp) >= 4) {
  sect("Fisher: adhesion vs proliferation x shared vs metastasis-only")
  print(table(mp$category, mp$sharing)); print(fisher.test(table(mp$category, mp$sharing)))
}

cat("\nNothing was written to disk. Every table is printed above and is also in\n",
    "these objects: cluster_proportions, cluster_genes, profile, excluded,\n",
    " gene_audit, matched_sets, model_counts, tree_meta, trees, mut_sites.\n",
    "Export whichever ones you want, e.g. write_csv(profile, \"profile.csv\")\n")
