# ------------------------------------------------------------------
# Masterfile 
# ------------------------------------------------------------------
rm(list = ls())

# Prerequisites:
# - These folders: inputs, outputs and R
# - All R files from github should be deposited in R folder
# - All packages should be installed from the setup.R script
# - Run the Individual_variant_analysis.R script once before the New_GRS_analysis.R script

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

# Age sensitivity analysis
source("R/Age_sensitivity_analyses.R")


# Go through scripts and make them tidier |+|
# Check that they all run in sequence |+|
# Make it so that setup is the only thing that needs to be ran and each of the analysis steps can be ran without running the other steps |+|
# Add in the age sensitivity analysis |+|
# Add in the metabolite GWAS into the age sensitivity analysis
# Change the age sensitivity analysis so it looks nice and clear
# Go through to see how genetic data and metabolite data is read in |+|
# Check the HMGCR step in setup - why does there need to be a variable HMGCR_SNPs?
# Maybe change the for loop in the lipoprotein groups bit - little confusing - could just label it better? |+|
# Change the heatmap production section so it produces all 4 in a smaller amount of code |+|
# Add new control check for biology-driven clustered heatmaps - at the moment uses "temp_dat" which is not good... |+|
# Find a way to add gene names to heatmaps without them going off the edge - then add gene names to all heatmaps
# Put linear regression function in setup.R |+|
# Check they all run in sequence again
# Remove prerequisites from Masterfile.R and put them into the README
# Need to change Dendrogram_production.R if keeping the different ages thing...





