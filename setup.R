##############
# data setup #
##############

#install.packages("cluster")
#install.packages("RColorBrewer")
#install.packages("gplots")
#install.packages("ape")
#install.packages("foreign")
#install.packages("openxlsx")
#install.packages("GenABEL")
#install.packages("gtools")
#install.packages("ggplot2")
#install.packages("tidyverse")
#install.packages("haven")
#install.packages("readxl")
#install.packages("qqman")
#install.packages("zeallot")
#install.packages("data.table")
#install.packages("gridExtra")

library(cluster)
library(RColorBrewer)
library(gplots)
library(ape)
library(foreign)
library(openxlsx)
library(GenABEL)
library(gtools)
library(ggplot2)
library(tidyverse)
library(haven)
library(readxl)
library(qqman)
library(stringr)
library(zeallot)
library(data.table)
library(gridExtra)

##################################################
### Read in all the data from Stata 	##########  #####NEED TO CHANGE THIS SECTION SO THAT YOU READ IN THE FILE PEOPLE WILL RECEIVE FROM ALSPAC AND CHANGE IT SO THAT IT IS MERGED - ALSO NEED TO CHANGE/ADD IN HOW TO SORT THE GENOTYPE FILE
##################################################

###Set working directory for file with everything in
#setwd("CAD_adolescent_analysis")

#the original data file with metabs
datafile_metabs <- read_dta("inputs/metabolite_data.dta") #was export_metabolomics

#the sorted SNP file
datafile_SNPs <- read_dta("inputs/genotype_data.dta") #was New_Genotype_Sorted_dosage_genotype_data
#HMGCR_SNPs <- read_dta("/Volumes/tb13101/Desktop/Mini-project 1/CAD_GRS_analysis/New_genotype_dosage_data/HMGCR/HMGCR_sorted_data.dta")
#datafile_SNPs <- datafile_SNPs %>%
 # mutate(u_ID = paste(cidB9999, qlet, sep = "_")) %>%
  #select(u_ID, everything())
#HMGCR_SNPs <- HMGCR_SNPs %>%
 # mutate(u_ID = paste(cidB9999, qlet, sep = "_")) %>%
  #select(u_ID, everything())

#write.dta(datafile_SNPs, "inputs/genotype_data.dta")

#the SNPs and their weightings
SNP_Ws <- read_dta("inputs/SNP_weightings.dta") #was SNPs and log_odd_weightings_new_genotype

SNP_info <- read_excel("inputs/SNP_info.xlsx")
colnames(SNP_info) <- c("Locus", "Lead_variant", "A1", "A2", "A1_freq", "OR", "logOR")
head(SNP_info)

#Add in to read in the alleles and then pull the reverse SNPs from differences between this and 
#SNP_info 
genotype_alleles <- read.table("inputs/alleles.txt", header = F, stringsAsFactors = F)
colnames(genotype_alleles) <- c("CHR", "SNP", "POS", "other_allele", "ref_allele")
#Change indels to insertion/deletion rather than actual bases
for (i in 1:nrow(genotype_alleles)) {
  if (nchar(genotype_alleles[i, "ref_allele"]) > 1) {
    genotype_alleles[i, "ref_allele"] <- "I"
    genotype_alleles[i, "other_allele"] <- "D"
  } else if (nchar(genotype_alleles[i, "other_allele"]) > 1) {
    genotype_alleles[i, "other_allele"] <- "I"
    genotype_alleles[i, "ref_allele"] <- "D"
  }
}

HMGCR_SNPs <- c("rs17238484", "rs12916")
#SNP_info <- filter(SNP_info, !(Lead_variant %in% HMGCR_SNPs))

SNP_info <- arrange(SNP_info, Lead_variant)
genotype_alleles <- arrange(genotype_alleles, SNP)
stopifnot(SNP_info$Lead_variant == genotype_alleles$SNP)

#Extract SNPs where the effect allele does not equal the reference allele
rev_SNPs <- genotype_alleles %>%
  arrange(SNP) %>%
  mutate(Correct_ref = ifelse(ref_allele == SNP_info$A1, "T", "F")) %>%
  filter(Correct_ref == "F") %>%
  .[["SNP"]]
rev_SNPs <- c(rev_SNPs, HMGCR_SNPs)

#########################################################
###  Produce SNP weightings  and sort out the SNPs  ### 
#########################################################

##NB: Weightings are based on the log(OR) of each SNP with regards to CAD

#produce a new SNP file that can be altered and remove the IDs
#datafile_SNPs_W <- datafile_SNPs[-c(1,2)]
datafile_SNPs_W <- datafile_SNPs
#str(datafile_SNPs_W)
#str(SNP_Ws)

#This loop changes the coding of the alleles so that the dosages all reflect the dosage of effect
#allele received by each individual then weights each SNP dosage
other_headers <- c("u_ID", "cidB9999", "qlet", HMGCR_SNPs)
headers <- colnames(datafile_SNPs_W)[!(colnames(datafile_SNPs_W) %in% other_headers)]
for (i in headers) {
  #Changes allele codings
  if (i %in% rev_SNPs) {
    datafile_SNPs_W[[i]] <- 2 - datafile_SNPs_W[[i]]
  }
  #Weights the SNPs
  weight <- subset(SNP_Ws, leadvariant == i)
  new_col_name <- paste(i, "_w", sep = "")
  datafile_SNPs_W[[new_col_name]] <- datafile_SNPs_W[[i]] * weight[[2]]
}
#datafile_SNPs_W[[HMGCR_SNPs[1]]] <- 2 - datafile_SNPs_W[[HMGCR_SNPs[1]]]
#datafile_SNPs_W[[HMGCR_SNPs[2]]] <- 2 - datafile_SNPs_W[[HMGCR_SNPs[2]]]
datafile_SNPs_W <- datafile_SNPs_W %>%
  mutate(u_ID = paste(cidB9999, qlet, sep = "_")) %>%
  dplyr::select(u_ID, everything())

#########################################################
## Merge the ages and data frames and sort out the dat ##
#########################################################
source("R/Merge.R")
#merge
d <- left_join(df_main, datafile_SNPs_W) %>%
  dplyr::select(u_ID, cidB9999, qlet, age, everything())
d

#extract SNP names - as well as lipid SNPs 
SNPs <-  colnames(d)[str_detect(colnames(d), "[\\d]_w$")] #try "\\d" in the square brackets

#save the name of the SNP that has a log(OR) of 0
log_odd_SNP <- "rs6903956_w"

#order the data by unique identifiers
d <- arrange(d, cidB9999, qlet)

length(mnames) - length(nr_mnames) #81

#For age specific analyses
if (age != "all") {
  d <- d[d[["age"]] == age, ]
  print(paste("The following analyses will be conducted at age", age, "only"))
}

if (age == 17) {
  d <- dplyr::select(d, -insulin, -glucose)
  nr_mnames <- nr_mnames[!(nr_mnames %in% c("insulin", "glucose"))]
}


print("setup complete - please proceed to the CAD-GRS analysis")




