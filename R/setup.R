# ------------------------------------------------------------------
# Data setup 
# ------------------------------------------------------------------

pkgs <- c("cluster", "RColorBrewer", "gplots", "ape", "foreign", "openxlsx", "GenABEL", "gtools", "ggplot2", "tidyverse", "haven", "readxl", "qqman", "zeallot", "data.table", "gridExtra", "FSA", "stringr")

for (i in pkgs) {
  if (!(i %in% installed.packages())) {
    install.packages(i)
  }
  library(i, character.only = T)
}

source("R/Linear_regression_func.R")

# For tables
make_pretty <- function (num, digits) {
  as.numeric(formatC(signif(num, digits), digits = digits, format = "fg", flag = "#"))
}

# ------------------------------------------------------------------
# Read in all the data 
# ------------------------------------------------------------------

datafile_metabs <- read_dta("inputs/metabolite_data.dta") 
#datafile_metabs <- datafile_metabs[, -grep(paste("glucose", "insulin", sep = "|"), colnames(datafile_metabs))]

datafile_SNPs <- read_dta("inputs/genotype_data.dta") # was New_Genotype_Sorted_dosage_genotype_data
datafile_SNPs <- datafile_SNPs %>%
  mutate(u_ID = paste(cidB9999, qlet, sep = "_")) %>%
  select(u_ID, everything())

SNP_info <- read_excel("inputs/SNP_info.xlsx")
colnames(SNP_info) <- c("Locus", "Lead_variant", "A1", "A2", "A1_freq", "OR", "logOR")
SNP_info

# ------------------------------------------------------------------
# Produce weighted variants with correct effect allele
# ------------------------------------------------------------------

genotype_alleles <- read.table("inputs/alleles.txt", header = F, stringsAsFactors = F)
colnames(genotype_alleles) <- c("CHR", "SNP", "POS", "other_allele", "ref_allele")

# Change indels to insertion/deletion rather than actual bases
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

SNP_info <- arrange(SNP_info, Lead_variant)
genotype_alleles <- arrange(genotype_alleles, SNP)
stopifnot(SNP_info$Lead_variant == genotype_alleles$SNP)

# Extract SNPs where the effect allele does not equal the reference allele
rev_SNPs <- genotype_alleles %>%
  arrange(SNP) %>%
  mutate(Correct_ref = ifelse(ref_allele == SNP_info$A1, "T", "F")) %>%
  filter(Correct_ref == "F") %>%
  .[["SNP"]]
rev_SNPs <- c(rev_SNPs, HMGCR_SNPs) 

## NB: Weightings are based on the log(OR) of each variant with regards to CAD

datafile_SNPs_W <- datafile_SNPs

# Change the coding of effect alleles where needed and weight the SNPs
other_headers <- c("u_ID", "cidB9999", "qlet", HMGCR_SNPs)
variants <- SNP_info$Lead_variant[!(SNP_info$Lead_variant %in% HMGCR_SNPs)]
for (i in variants) {
  # Changes allele codings
  if (i %in% rev_SNPs) {
    datafile_SNPs_W[[i]] <- 2 - datafile_SNPs_W[[i]]
  }
  # Weights the SNPs
  weight <- subset(SNP_info, Lead_variant == i)
  new_col_name <- paste(i, "_w", sep = "")
  datafile_SNPs_W[[new_col_name]] <- datafile_SNPs_W[[i]] * weight[["logOR"]]
}
datafile_SNPs_W <- datafile_SNPs_W %>%
  mutate(u_ID = paste(cidB9999, qlet, sep = "_")) %>%
  dplyr::select(u_ID, everything())

# ------------------------------------------------------------------
# Merge the ages and and extract useful info
# ------------------------------------------------------------------
source("R/Merge.R")
# merge 
d <- left_join(df_main, datafile_SNPs_W) %>%
  dplyr::select(u_ID, cidB9999, qlet, age, everything())
d

# Extract SNP names - as well as lipid SNPs 
SNPs <-  colnames(d)[str_detect(colnames(d), "[\\d]_w$")] 

# Removing the SNPs that didn't meet genome-wide significance in Nikpay et al. paper
#SNP_to_rm <- c("rs273909", "rs6903956", "rs17609940", "rs10953541", "rs264", "rs2954029", "rs964184", "rs9319428", "rs17514846", "rs216172", "rs12936587", "rs46522")
#SNP_to_rm <- paste0(SNP_to_rm, "_w")

#SNPs <- SNPs[!(SNPs %in% SNP_to_rm)]

#d <- d[, !(colnames(d) %in% SNP_to_rm)]



lipoproteins <- c(nr_mnames[grep("-V?[HIL]DL", nr_mnames)], nr_mnames[grep("IDL", nr_mnames)])
non_lipo <- nr_mnames[!(nr_mnames %in% lipoproteins)]

# Re-order the data by unique identifiers
d <- arrange(d, cidB9999, qlet)

length(mnames) - length(nr_mnames) #81

# For age specific analyses
if (age != "all") {
  d <- d[d[["age"]] == age, ]
  print(paste("The following analyses will be conducted at age", age, "only"))
}

if (age == 17) {
  d <- dplyr::select(d, -insulin, -glucose)
  nr_mnames <- nr_mnames[!(nr_mnames %in% c("insulin", "glucose"))]
}

# ------------------------------------------------------------------
# Generate biological subsets
# ------------------------------------------------------------------
# Generate the lipoprotein subsets 
Large_VLDL <- grep(paste(c("xxl-VLDL", "xl-VLDL", "l-VLDL"), collapse = "|"), nr_mnames, value = TRUE)
Remnant_particles <- grep(paste(c("m-VLDL", "s-VLDL", "xs-VLDL", "IDL"), collapse = "|"), nr_mnames, value = TRUE)
LDL <- grep(paste(c("l-LDL", "m-LDL", "s-LDL"), collapse = "|"), nr_mnames, value = TRUE)
V_Large_HDL <- grep(c("xl-HDL"), nr_mnames, value = TRUE) 
Large_HDL <- grep(paste(c("^l-HDL", "m-HDL"), collapse = "|"), nr_mnames, value = TRUE)
Small_HDL <- grep("s-HDL", nr_mnames, value = TRUE)

# Generate the other metabolite subsets
diameter <- non_lipo[grep("-D", non_lipo)]
cholesterol <- non_lipo[grep("-C$", non_lipo)]
glycerides <- non_lipo[grep("-TG$", non_lipo)]
phospholipids <- c("DAG", "Tot-PG", "PC", "Tot-Cho")
apolipo <- non_lipo[grep("Apo", non_lipo)]
fatty_acids <- c("Tot-FA", "FALen", "UnsatDeg", "DHA", "LA", "CLA", "FAw3", "FAw6", "PUFA", "MUFA", "SFA")
amino_acids <- c("Ala", "Gln", "His", "Ile", "Leu", "Val", "Phe", "Tyr")
glycolysis <- c("Glc", "Lac", "Pyr", "Cit", "glucose", "insulin")
other <- c("Ace", "AcAce", "bOHBut", "Crea", "Alb", "Gp")

subsets <- list(Large_VLDL, Remnant_particles, LDL, V_Large_HDL, Large_HDL, Small_HDL, diameter, cholesterol, glycerides, phospholipids, apolipo, fatty_acids, amino_acids, glycolysis, other)
names(subsets) <- c("Large_VLDL", "Remnant_particles", "LDL", "V_Large_HDL", "Large_HDL", "Small_HDL", "diameter", "cholesterol", "glycerides", "phospholipids", "apolipo", "fatty_acids", "amino_acids", "glycolysis", "other")

subset_df <- data.frame(Metabolite = nr_mnames, subset = NA, group = NA)

for (i in names(subsets)) {
  subset_df[subset_df[["Metabolite"]] %in% get(i), "subset"] <- i
}

# Generate the lipoprotein groups
est_tot <- lipoproteins[grep(pattern = "-C", lipoproteins)]
free <- lipoproteins[grep(pattern = "-FC", lipoproteins)]
cholesterol <- lipoproteins[lipoproteins %in% c(free, est_tot)]
phospholipid <- lipoproteins[grep(pattern = "-PL", lipoproteins)]
triglyceride <- lipoproteins[grep(pattern = "-TG", lipoproteins)]
tot_lipid <- lipoproteins[grep(pattern = "L-L", lipoproteins)]
particle <- lipoproteins[grep(pattern = "-P", lipoproteins)]
particle <- particle[-which(particle %in% phospholipid)]

groups <- list(cholesterol, phospholipid, triglyceride, tot_lipid, particle)
names(groups) <- c("cholesterol", "phospholipid", "triglyceride", "tot_lipid", "particle")

for (i in names(groups)) {
  subset_df[subset_df[["Metabolite"]] %in% groups[[i]], "group"] <- i
}

write.table(subset_df, file = "outputs/tables/subset_table.txt", quote = F, col.names = T, row.names = F, sep = "\t")

print("setup complete - please proceed to the CAD-GRS analysis")




