# ------------------------------------------------------------------
# Masterfile 
# ------------------------------------------------------------------
rm(list = ls())

# ------------------------------------------------------------------
# Setup 
# ------------------------------------------------------------------

# source("R/setup.R")
source("R/new_setup.R")

# ------------------------------------------------------------------
# Analysis 
# ------------------------------------------------------------------

# Analysis of the association between CAD-GRS and metabolites
source("R/CAD_GRS_analysis.R") 

# Analysis of the association between all variants and metabolites
source("R/Individual_variant_analysis.R")

# Analysis of association between SNPs within HMGCR (proxies for statin use) and metabolites
source("R/HMGCR_SNP_analysis.R")

# Age sensitivity analysis - STILL WORKING ON IT! 
source("R/Age_sensitivity_analyses.R")





