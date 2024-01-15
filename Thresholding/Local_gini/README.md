# LOCAL GINI THRESHOLDING

Local T2 thresholding method relies on local criteria to determine whether a gene is "core" or "non-core" in a given sample. LocalT2 adds upper and lower bounds to the defining of thresholds, making it easier to distinguish between core and non-core genes.
## Components
The folder includes different files.

- **Data**: Merged_data, Mod_data and NM2ENSG.

- **Metabolic models**: Human-GEM_Cobra_v1.01, SysBio_COBRA_v1.13 and SysBio_COBRA_v1.17_consensus.
- **Functions**: GiniReactionImportance.m buildContextmodels.m.
- **Results**: The code takes a lot of time to run, so I already uploaded the results too (if you want to run by yourself, check the code, and comment the lines that load it, and remove the coments of those that call the function), RxnImp.mat, reactions_hkg_hkg.mat

- **Localgini_threshold**: This code sets the core genes based on the local gini definition done by Kumar & Bhatt., 2023, explained below. Once obtained the core genes, they are mapped to the reactions they are involved, giving as result the core reactions. The core genes and core reactions are compared with housekeeping genes and reactions respectively.
## Local gini theory
Localgini is a bioinformatics technical method that uses the Gini coefficient to determine thresholds in gene expression analysis. Localgini, in contrast to classic approaches, assigns specific thresholds for each gene based on its expression distribution across different samples. This strategy enables for improved capture of fluctuations in gene expression, allowing for the inclusion of key genes in genomic models that would otherwise be excluded by other methods. Because of its emphasis on gene-specific thresholds, it is particularly effective for developing more accurate and representative genomic models.
## Bibliography
The Local gini thresholding code was written following the considerations and the github repository of Kumar & Bhatt, 2023.
- Kumar S, P., & Bhatt, N. (2023). Localgini: A method for harnessing inequality in gene expression to improve the quality of context-specific models. bioRxiv, 2023-09.
## Authors

- [@pablots98](https://www.github.com/pablots98)
