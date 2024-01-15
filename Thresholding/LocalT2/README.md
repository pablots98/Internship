
# LOCALT2 THRESHOLDING

Local T2 thresholding method relies on local criteria to determine whether a gene is "core" or "non-core" in a given sample. LocalT2 adds upper and lower bounds to the defining of thresholds, making it easier to distinguish between core and non-core genes.
## 
The folder includes different files.

- **Data**: Merged_data, Mod_data and NM2ENSG.

- **Metabolic models**: Human-GEM_Cobra_v1.01, SysBio_COBRA_v1.13 and SysBio_COBRA_v1.17_consensus.
- **Functions**: localT2_new.

- **Local_T2**: This code sets the core genes based on the local T2 definition done by Richelle et al., 2019, explained below. Once obtained the core genes, they are mapped to the reactions they are involved, giving as result the core reactions. The core genes and core reactions are compared with housekeeping genes and reactions respectively.
## LocalT2 theory
The localT2 processes gene expression data using lower and upper percentile thresholds. It converts the data to a binary format based on these thresholds, assigning 1 to expression values that meet or exceed the threshold and 0 to those that do not. The function adjusts the thresholds for each gene individually, based on its mean expression, and generates two matrices: one binary and one with adjusted expression values between 0 and 1. This is useful for normalising gene expression data for further analysis.
## 
The Local_T2 thresholding code was written following the considerations of Richelle et al., 2019. 
- Richelle, A., Joshi, C., & Lewis, N. E. (2019). Assessing key decisions for transcriptomic data integration in biochemical networks. PLoS computational biology, 15(7), e1007185.
##
## Authors

- [@pablots98](https://www.github.com/pablots98)
