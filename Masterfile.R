# ------------------------------------------------------------------
# Masterfile 
# ------------------------------------------------------------------
rm(list = ls())

# Prerequisites:
# - These folders: inputs, outputs and R
# - All R files from github should be deposited in R folder
# - All packages should be installed from the setup.R script



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

# Analysis of association between CAD-GRS and metabolites
source("R/CAD_GRS_analysis.R")

# Analysis of association between all variants and metabolites
source("R/Individual_variant_analysis.R")

# Analysis of association between new GRS (formed from SNPs that associate with metabolites) and metabolites
# with comparison to original CAD-GRS
source("R/New_GRS_analysis.R")

# Analysis of association between SNPs within HMGCR (proxies for statin use) and metabolites
source("R/HMGCR_SNP_analysis.R")


# Go through scripts and make them tidier
# Check that they all run in sequence
# Make it so that setup is the only thing that needs to be ran and each of the analysis steps can be ran without running the other steps
# Add in the age sensitivity analysis
# Add in the metabolite GWAS into the age sensitivity analysis
# Go through to see how genetic data and metabolite data is read in
# Don't want metabolite data that is currently read in to be what Hash did?
# Check the HMGCR step in setup - why does there need to be a variable HMGCR_SNPs?










