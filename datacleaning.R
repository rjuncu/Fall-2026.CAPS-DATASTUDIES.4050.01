#install.packages("readr")
#install.packages("tidyverse")

library(readr)
library(tidyverse)

clinical_data <- read_tsv("./coadread_mskcc_clinical_data.tsv")
colnames(clinical_data)

#Check that all Sample IDs are unique
clinical_data %>%
  count(`Sample ID`) %>%
  filter(n > 1)

#Extract the MSK cBioPortal archive
untar("coadread_mskcc.tar.gz")
list.files()
mutation_data <- read_tsv("./coadread_mskcc/data_mutations.txt")
colnames(mutation_data)

#Make sure clinical sampleID and mutation tumor_sample_barcode are the same dataset
#this will be the link between the two data sets 
head(unique(clinical_data$`Sample ID`))
length(unique(clinical_data$`Sample ID`))

head(unique(mutation_data$Tumor_Sample_Barcode))
length(unique(mutation_data$Tumor_Sample_Barcode))

common_samples <- intersect(
  clinical_data$`Sample ID`,
  mutation_data$Tumor_Sample_Barcode
)

length(common_samples)
#There are 138 unique samples shared between the two datasets

#Now 138 samples must be sorted into 69 patient pairs 
#Pairs include a primary tumor and a metastatic tumor
patient_pairs <- clinical_data %>%
  select("Patient ID", "Sample ID", "Sample Type") %>%
  tidyr::pivot_wider(
    names_from = "Sample Type",
    values_from = "Sample ID"
  )
head(patient_pairs)  

#check that every patient has two tumor samples
table(
  !is.na(patient_pairs$Primary),
  !is.na(patient_pairs$Metastasis)
)


#integrate mutation and clinical data by matching barcodes to sampleID
annotated_mutations <- mutation_data %>%
  filter(Tumor_Sample_Barcode %in% common_samples) %>%
  left_join(
    clinical_data %>%
      select("Patient ID", "Sample ID", "Sample Type"),
    by = c("Tumor_Sample_Barcode" = "Sample ID")
  )

colnames((annotated_mutations))

#Remove columns with no data observations 
annotated_mutations <- annotated_mutations %>%
  select(-Entrez_Gene_Id, -dbSNP_RS, -dbSNP_Val_Status, -Match_Norm_Seq_Allele1, -Match_Norm_Seq_Allele2, -Tumor_Validation_Allele1, -Tumor_Validation_Allele2, -Match_Norm_Validation_Allele1, -Match_Norm_Validation_Allele2, -Verification_Status, -Sequencing_Phase, -Sequence_Source, -Validation_Method, -Score, -BAM_File, -Sequencer)

#Remove columns with only one unique response 
#This will reduce dimensionality 
#The removed data is qualitative and provides context for the analysis 
#This context will be explored in the project background
annotated_mutations <- annotated_mutations %>%
  select(-Center, -NCBI_Build, -Strand, -Matched_Norm_Sample_Barcode, -Validation_Status, -Mutation_Status)

colnames(annotated_mutations)

#Check that there is no missing data 
annotated_mutations %>% summarise(
  total_mutations = n(), 
  unique_mutations = n_distinct(Tumor_Sample_Barcode), 
  unique_patients = n_distinct(`Patient ID`), 
  missing_patient_ids = sum(is.na(`Patient ID`)),
  missing_sample_types = sum(is.na(`Sample Type`))
)

#Calculate VAF
annotated_mutations <- annotated_mutations %>%
  mutate(VAF = t_alt_count/(t_ref_count + t_alt_count))

#Check that values are between 0 and 1
summary(annotated_mutations$VAF)

#Export into file for SciClone Analysis 
write_tsv(annotated_mutations, "annotated_mutations.tsv")

list.files()
