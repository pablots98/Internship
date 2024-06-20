# Script:       Step2_EnrichmentAnalysis.R
# Description:  In this script, we will explore the differential gene 
#               expression dataset comparing lung cancer vs. healthy tissue
#               samples. The RNA-sequencing dataset was retrieved from 
#               TCGA (The Cancer Genome Atlas) and pre-processed in R. 
#               Differential gene expression analysis was performed with the 
#               DESeq2 R-package. 
# Version: 1.0
# Last updated: 2024-06-10
# Author: mkutmon

# #############################################
# R INSTRUCTIONS
# #############################################

# Make sure you ran Step1_DataExploration.R directly before starting this script
library(clusterProfiler)
library(enrichplot)
# #############################################
# GENE ONTOLOGY ENRICHMENT ANALYSIS
# #############################################

# Load the data
Eis_met = read.table('unique_eis_met.txt', header = FALSE, sep = "\t")
HK_met = read.table('unique_hk_met.txt', header = FALSE, sep = "\t")
Joshi_met = read.table('unique_joshi_met.txt', header = FALSE, sep = "\t")

# We will first explore what kind of biological processes are affected by 
# performing a Gene Ontology enrichment analysis. 

out.folder = "C:/Users/PC/OneDrive/Documentos/Systems_Biology_master/Internship/Internship/Thresholding/GO_analysis"

# We will start by looking at processes that are up-regulated
res.go.eis <- clusterProfiler::enrichGO(Eis_met$V1, OrgDb = "org.Hs.eg.db", 
                   keyType="ENSEMBL", ont = "BP", 
                   pAdjustMethod = "BH", qvalueCutoff = 0.05, 
                   minGSSize = 20, maxGSSize = 400)
res.go.eis.df <- as.data.frame(res.go.eis)
write.table(res.go.eis.df, file=paste0(out.folder,"go-eis.txt"), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

res.go.eis.sim <- enrichplot::pairwise_termsim(res.go.eis)

# then we visualize them in a treeplot
treeplot(res.go.eis.sim, label_format = 0.5, showCategory = 100, cluster.params = list(n = 10))

# we will also save files that has a nicely readable figure
filename <- paste0(out.folder,"GO-treeplot-eis.png")
png(filename , width = 3000, height = 4000, res = 150)
plot(treeplot(res.go.eis.sim, label_format = 0.5, showCategory = 100, cluster.params = list(n = 10)))
dev.off()
