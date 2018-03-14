# ------------------------------------------------------------------
# Masterfile 
# ------------------------------------------------------------------
rm(list = ls())

# Select the age at which the analysis is to be conducted at
# If all ages make age "all"
age <- "all"
#age <- 7
#age <- 15
#age <- 17

# ------------------------------------------------------------------
# Setup 
# ------------------------------------------------------------------

source("R/setup.R")

# ------------------------------------------------------------------
# Analysis 
# ------------------------------------------------------------------

# Analysis of the association between CAD-GRS and metabolites
source("R/CAD_GRS_analysis.R") 

# Analysis of the association between all variants and metabolites
source("R/Individual_variant_analysis.R")

# Analysis of the association between the new GRS (formed from SNPs that associate with metabolites)
# and metabolites with comparison to original CAD-GRS
source("R/New_GRS_analysis.R")

# Analysis of association between SNPs within HMGCR (proxies for statin use) and metabolites
source("R/HMGCR_SNP_analysis.R")

# Age sensitivity analysis - STILL WORKING ON IT! 
source("R/Age_sensitivity_analyses.R")


# Add in the metabolite GWAS into the age sensitivity analysis
# Change the age sensitivity analysis so it looks nice and clear
# Find a way to add gene names to heatmaps without them going off the edge - then add gene names to all heatmaps
# Check they all run in sequence again 
# Need to change Dendrogram_production.R if keeping the different ages thing...





