# Fall-2026.CAPS-DATASTUDIES.4050.01
## Primary and Metastatic Colorectal Adenocarcinoma Tumor Reconstruction and Evolution
### Using SciClone and ClonEvol 
#### Data Source: https://www.cbioportal.org/study/clinicalData?id=coadread_mskcc 
Brannon AR, Vakiani E, Sylvester BE, Scott SN, McDermott G, Shah RH, Kania K, Viale A, Oschwald DM, Vacic V, Emde AK, Cercek A, Yaeger R, Kemeny NE, Saltz LB, Shia J, D'Angelica MI, Weiser MR, Solit DB, Berger MF. Comparative sequencing analysis reveals high genomic concordance between matched primary and metastatic colorectal cancer lesions. Genome Biol. 2014 Aug 28;15(8):454. doi: 10.1186/s13059-014-0454-7. PMID: 25164765; PMCID: PMC4189196.

#### **Files:** 

```coadread_mskcc.tar.gz```
> A compressed file containing mutational data matrices 

```coadread_mskcc_clinical_data.tsv```
> Clinical data that was scraped and compressed from ```coadread_mskcc.tar.gz```

```annotated_mutations.tsv```
> Merged clinical and mutational data that has been prepped for input into SciClone

```sciClone_testcase.R```
> Script creating subclonal clusters from sequencing read data in ```annotated_mutations.tsv```

```sciClone_full_analysis.R```
> Script for applying SciClone clustering on pairs of patient segments that were filtered in ```sciClone_testcase.R```

```sciclone_patient_input.tsv```
> Dataset containing a list of patients who have at least 3 shared mutations between primary and metastatic tumors, with a depth of at least 100

```sciclone_mutation_input.tsv```
> Dataset containing mutation data from patients listed in ```sciclone_patient_input.tsv```
