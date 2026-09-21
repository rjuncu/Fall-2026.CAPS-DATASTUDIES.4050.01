# Fall-2026.CAPS-DATASTUDIES.4050.01
## Primary and Metastatic Colorectal Adenocarcinoma Tumor Reconstruction and Evolution
### Using SciClone and ClonEvol 

**Files:** 

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
