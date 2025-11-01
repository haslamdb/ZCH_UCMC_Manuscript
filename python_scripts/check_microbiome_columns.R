#!/usr/bin/env Rscript
# Quick script to check what columns are in NICUSpeciesNR

load("/home/david/projects/Metagenomics/Yanping/NICU_Microbiome/Hangzhou/NoHumanDNA20220929/NICUData20250514")

cat("Column names in NICUSpeciesNR:\n")
print(head(colnames(NICUSpeciesNR), 20))

cat("\n\nFirst few rows of key columns:\n")
print(head(NICUSpeciesNR[, 1:10]))

cat("\n\nUnique values in first column:\n")
print(head(unique(NICUSpeciesNR[, 1]), 20))
