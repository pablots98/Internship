# Internship
## Context-Specific Model for Senescent Fibroblasts
This repository contains the development of a specific contextual model to study fibroblasts in senescence. It aims to analyse changes in fibroblasts during ageing and their impact on tissues and organs. It includes code, data and related documentation.
## 
- **Data**: Includes all the data, from the transcriptomics dataset to the list of housekeeping genes

- **FPKM to TPM**: The data from the dataset was normalized in FPKM. Based on different bibliographic research done by Shuyi, discovered that TPM actually can be a better way of normalization to the comparison across samples.

- **Thresholding**: To run mCADRE or other algorithms to create Context-Specific Models, it is necessary to create a vector or matrix with the core reactions based on the transcriptomics data. Here we study different thresholding methods to establish which one performs better. To investigate that, we use the housekeeping genes and reactions, based on the assumption that the housekeeping reactions should be considered as core reactions. Based on that assumption, the thresholding method that includes more housekeeping reactions as core reactions is going to be the most accurate one.

- **mCADRE**:
