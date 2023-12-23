# Gene-Expression-Analysis
Bio Principle Final Assignment
# Usage
Download the file ，open final.Rproj, and run the **Gene Expression Analysis.r**


The data would be untared and set workpath to the data file.
s
# Different Part of Code

## Import Data, Preprocess, Create Metadata

1. Download the dataset on: [cbioportal](https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018)

2.  Untar the folder and extract the
   files.

3.  Read the RNASeq file:
   data_mrna_seq_v2_rsem.txt

4. Read the Patient Data file:
   data_clinical_patient.txt

5. Read the Copy Number Aberrations
   Data: data_cna.txt

6. Match the RNASeq patient ids with
   the CNA ids and the Patient Data ids.

7. Create metadata using the CNA level of ERBB2+
   (greater than 0 means amplified).

## 1. Differential Expression Analysis

**HER2 Amplified and Not Amplified**

Normalize data using DESeq2.

Obtain Differentially Expressed Genes.

## 2. Top 10 Differentially Expressed Genes

**Ranked by Fold Change**

## 3. Pathway Enrichment

## 4. Get The Variance Stabilized Transformed Expression Values

## 5.PCA Plot

With the vst values obtain a PCA plot.

## 6. Gene Expression Clusters

Cluster the data and show in PCA

## 7. Cox Survival Regression model

### 1.   With the vst values of the DE genes generate an overall survival model

### 2. Use lasso cross validation to find the set of genes which predict survival the best.






