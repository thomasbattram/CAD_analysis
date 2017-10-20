######################
# Masterfile #########
######################
rm(list = ls())

#########
# Setup #
#########
source("R/setup.R")

############
# Analysis #
############

#Analysis of association between CAD-GRS and metabolites
source("R/CAD_GRS_analysis.R")

#Analysis of association between all variants and metabolites
source("R/Individual_variant_analysis.R")

#Analysis of association between new GRS (formed from SNPs that associate with metabolites) and metabolites
#with comparison to original CAD-GRS
source("R/New_GRS_analysis.R")

#Analysis of association between SNPs within HMGCR (proxies for statin use) and metabolites
source("R/HMGCR_SNP_analysis.R")


#All ggplot graphs - larger writing
#Heatmaps - add in gene names
