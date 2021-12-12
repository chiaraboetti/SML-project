# Gene selection for cancer classification
Final project of SML course a.y. 2021-2022 at Stochastic and Data Science (UniTo).  Project due to Friday, December 17 (6pm).

## The problem
The task is Gene Selection for Cancer Classification, i.e. finding a reasonably low number of important genes and use them to classify the type of cancer. To this aim, three learning algorithms were used:
- Random Forests
- SVM-Lasso
- Feed Forward Neural Networks

In SVM-Lasso, important features are considered to be the ones which have a non-zero weight. In the other two cases, important genes were selected after fitting the global model and using spefic procedures of Variable Importance such as Permutation Imprtance and Olden Importance.

## Datasets
Data files can be found on the Depmap website.  This project aims to create  “cancer dependency map” by systematically identifying genetic dependencies and small molecule sensitivities and discovering the biomarkers that predict them.
In particular two datasets are used:

- CRISPR_gene_dependency.csv, which contains a sample of 1032 cancer cell lines, each of which is characterised by 17393 consecutive genes.

- sample_info.csv, which contains the cell line information (in particular the type of cancer).

The latter can be found in the dataset folder, the former is too big to be uploaded but can be downloaded freely from https://depmap.org/portal/download/ (looking for DepMap Public 21Q3, released on August, 2021).
