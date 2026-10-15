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

sciclone_total <- read_tsv("./SciClone Analysis/sciclone_results/tables/sciclone_total")
summary(sciclone_total)

ev001 <- sciclone_total %>%
  filter(
    Patient_ID == "EV-001",
    adequateDepth == 1,
    !is.na(cluster),
    cluster > 0,
    !is.na(Primary_VAF),
    !is.na(Metastasis_VAF)
  ) %>%
  select(
    cluster, 
    Primary = Primary_VAF,
    Metastasis = Metastasis_VAF
  ) %>%
  arrange(cluster)

print(ev001)
summary(ev001)

library(dplyr)

# Classify each observation without changing its original assignment
diagnostic <- sciclone_total %>%
  mutate(
    assignment_status = case_when(
      is.na(adequateDepth) ~ "Unknown depth status",
      adequateDepth != 1 ~ "Inadequate depth",
      is.na(cluster) ~ "Adequate depth, cluster NA",
      cluster == 0 ~ "Adequate depth, cluster 0",
      cluster > 0 ~ "Adequate depth, assigned cluster",
      TRUE ~ "Other"
    )
  )

# Overall counts
diagnostic %>%
  count(assignment_status, sort = TRUE)

# Compare assignment status across patients
diagnostic %>%
  count(Patient_ID, assignment_status) %>%
  arrange(Patient_ID, assignment_status)

# Check how many positive clusters each patient has
diagnostic %>%
  filter(
    adequateDepth == 1,
    !is.na(cluster),
    cluster > 0
  ) %>%
  group_by(Patient_ID) %>%
  summarise(
    n_mutations = n(),
    n_clusters = n_distinct(cluster),
    cluster_ids = paste(sort(unique(cluster)), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(n_clusters, Patient_ID)


# Count cluster-0 observations for each patient
cluster0_by_patient <- diagnostic %>%
  filter(adequateDepth == 1, cluster == 0) %>%
  group_by(Patient_ID) %>%
  summarise(
    n_cluster0 = n(),
    mean_primary_vaf = mean(Primary_VAF, na.rm = TRUE),
    mean_metastasis_vaf = mean(Metastasis_VAF, na.rm = TRUE),
    min_primary_depth = min(Primary_depth, na.rm = TRUE),
    min_metastasis_depth = min(Metastasis_depth, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_cluster0))

print(cluster0_by_patient, n = Inf)


diagnostic %>%
  filter(adequateDepth == 1) %>%
  group_by(Patient_ID, cluster) %>%
  summarise(
    n_mutations = n(),
    primary_vaf_mean = mean(Primary_VAF, na.rm = TRUE),
    metastasis_vaf_mean = mean(Metastasis_VAF, na.rm = TRUE),
    primary_vaf_sd = sd(Primary_VAF, na.rm = TRUE),
    metastasis_vaf_sd = sd(Metastasis_VAF, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Patient_ID, cluster)

# Identify mutations using patient + chromosome + position
cluster0_check <- diagnostic %>%
  filter(adequateDepth == 1, cluster == 0) %>%
  mutate(mutation_key = paste(Patient_ID, chr, st, sep = "_"))

clustered_check <- sciclone_total %>%
  mutate(mutation_key = paste(Patient_ID, chr, st, sep = "_"))

# Which cluster-0 mutations are absent from the clustered table?
cluster0_missing_from_clean <- cluster0_check %>%
  anti_join(
    clustered_check %>% select(mutation_key),
    by = "mutation_key"
  )

nrow(cluster0_missing_from_clean)
print(cluster0_missing_from_clean, n = Inf)

# shorten vaf column names as they will be
vaf.col.names <- grep('.vaf', colnames(x), value=T)
sample.names <- gsub('.vaf', '', vaf.col.names)
x[, sample.names] <- x[, vaf.col.names]
vaf.col.names <- sample.names
# prepare sample grouping
sample.groups <- c('P', 'R');
names(sample.groups) <- vaf.col.names
# setup the order of clusters to display in various plots (later)
x <- x[order(x$cluster),]