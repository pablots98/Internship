# StanDep THRESHOLDING

It employs a statistical method to evaluate gene reliance on various situations, allowing for more accurate modeling of biochemical networks. StanDep focuses on identifying genes whose expression is required for specific metabolic activities and distinguishes itself by its capacity to effectively integrate huge gene expression datasets, hence leading to a deeper understanding of biological systems.
## Components
The folder includes different files.

- **Data**: Merged_data, Mod_data and NM2ENSG.

- **Metabolic models**: Human-GEM_Cobra_v1.01, SysBio_COBRA_v1.13 and SysBio_COBRA_v1.17_consensus.
- **Functions**: clusterVariability1.m comparePromiscuousSpecific.m geneExprDist_hierarchy.m getModelData.m getPromEnzymes.m getUbiquityScore_2022, linearization_index.m models4mClusters1.m.

- **Project_StanDep**: This code sets the core reactions based on the StanDep definition done by Joshi et al., 2020, explained below. The core reactions are compared with housekeeping genes and reactions respectively.
## StanDep theory
StanDep is a sophisticated method for incorporating transcriptomic data into metabolic models at the genomic level. It employs a novel approach that employs certain thresholds for each group of genes, based on the variability of their expression under various conditions. This strategy allows for the capture of a broader range of genes related to cell function, including those with low expression yet important biologically. StanDep significantly improves the quality and precision of contextual metabólic models by taking into account the complexity and variability of genetic expression.
## 
The StanDep thresholding code was written following the considerations and the github repository of Joshi et al., 2020.
- Joshi, C. J., Schinn, S. M., Richelle, A., Shamie, I., O’Rourke, E. J., & Lewis, N. E. (2020). StanDep: Capturing transcriptomic variability improves context-specific metabolic models. PLoS computational biology, 16(5), e1007764..
##
## Authors

- [@pablots98](https://www.github.com/pablots98)
