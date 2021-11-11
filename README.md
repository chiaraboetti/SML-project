# Gene selection for cancer classification
Final project of SML course a.y. 2021-2022 at Stochastic and Data Science (UniTo).  Project due to Friday, December 17 (6pm) and review received on Monday, December 20 and due by Tuesday, December 21 (midnight).
The goal is to solve a real problem with a family of algorithms studied in class (or a new one).

## The problem
In this work, we applied machine learning algorithms to find the optimal set to classify cancer types with raw data. The hypothesis is that each cancer type has a set of genes representing the whole growth phenotype of such specific cancer.
Our primary goal is to identify a minimal feature (gene) set that still achieves reasonable classification of cancer type, so that feature selection is the key problem in this project.

## Datasets
Data files can be found on the Depmap website.  This project aims to create  “cancer dependency map” by systematically identifying genetic dependencies and small molecule sensitivities and discovering the biomarkers that predict them.
In particular two datasets are used:

- CRISPR_gene_dependency.csv, which contains a sample of 1032 cancer cell lines, each of which is characterised by 17393 consecutive genes.

- sample_info.csv, which contains the cell line information (in particular the type of cancer).

The latter can be found in the dataset folder, the former is too big to be uploaded but can be downloaded freely from https://depmap.org/portal/download/ (looking for DepMap Public 21Q3, released on August, 2021).
