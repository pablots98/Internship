setwd("C:/Users/PC/OneDrive/Documentos/Systems_Biology_master/Internship/Internship/Thresholding/HKgenes")

# Instalar y cargar biomaRt
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

library(biomaRt)

# Conectar con la base de datos Ensembl

# Lista de ENSEMBL Transcript IDs
data <- read.csv("Housekeeping_GenesHuman.csv", sep = ";")


transcript_ids <- data$Ensembl


mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")


res <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'ensembl_gene_id', 
                            'external_transcript_name',
                            'external_gene_name'),
             filters = 'ensembl_transcript_id_version', 
             values = transcript_ids,
             mart = mart)
