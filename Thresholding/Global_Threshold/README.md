GLOBAL THRESHOLDING
Global thresholding is the simplest thresholding method. It works by analysing all the genes in each sample, based on their data it establishes a percentile (you can set the one you want), where the genes that are expressed above will be considered core-genes, and those that are not, considered non-core.

The folder includes different files.

Data: Merged_data, Mod_data and NM2ENSG.

Metabolic models: Human-GEM_Cobra_v1.01, SysBio_COBRA_v1.13 and SysBio_COBRA_v1.17_consensus.

Global_Threshold_reactions: This code sets the core genes based on the above, and studies the percentage of core genes that are considered housekeeping genes. In turn, through these core-genes, the core-reactions needed to run the mCADRE algorithm are established.

The Global thresholding code was written following the considerations of Richelle et al., 2019.

Richelle, A., Joshi, C., & Lewis, N. E. (2019). Assessing key decisions for transcriptomic data integration in biochemical networks. PLoS computational biology, 15(7), e1007185.
Authors
@pablots98
