################### Biomart Thing ##############################################
################################################################################

# Clean the environment
rm(list = ls())
cat("\014")

## set the libraries
library(biomaRt)

# Set the working directory
setwd("C:/Users/PC/OneDrive/Documentos/Systems_Biology_master/Internship/Internship/Data")

## Load the data 
h_k_g = read.csv('Housekeeping_GenesHuman.csv', header = TRUE, sep = ';', stringsAsFactors = FALSE)

# Separate the ENS T ID
h_k_g_ENST = h_k_g$Ensembl
print(h_k_g_ENST)



################################################################################
##                            BIOMART THING                                   ##
################################################################################
