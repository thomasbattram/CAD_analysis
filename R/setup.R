# ------------------------------------------------------------------
# Data setup 
# ------------------------------------------------------------------

pkgs <- c("cluster", "RColorBrewer", "gplots", "ape",
          "foreign", "openxlsx", "GenABEL", "gtools",
          "ggplot2", "tidyverse", "haven", "readxl", 
          "qqman", "zeallot", "data.table", "gridExtra",
          "FSA", "stringr", "cowplot", "ggdendro")
lapply(pkgs, require, character.only = T)

source("R/Linear_regression_func.R")

# For tables
make_pretty <- function (num, digits) {
  as.numeric(formatC(signif(num, digits), digits = digits, format = "fg", flag = "#"))
}

args <- commandArgs(trailingOnly = T)
met_file <- args[1]
gen_file <- args[2]
cohort_char_file <- args[3]

# ------------------------------------------------------------------
# Read in all the data 
# ------------------------------------------------------------------
datafile_metabs <- read_dta(met_file) %>%
  mutate(u_ID = paste(aln, qlet, sep = "_")) %>%
  dplyr::select(u_ID, everything())

# Set missing values to NA and withdrawn consent to NA
datafile_metabs[,-c(1:4)] <- apply(datafile_metabs[,-c(1:4)], 2, function(x) {replace(x, which(x == -1 | x == -9999), NA)})

# remove pyruvate because it is indistinguishable from EDTA in NMR analysis
datafile_metabs <- datafile_metabs[, -grep("pyr", colnames(datafile_metabs), ignore.case = T)]

datafile_SNPs <- read.table(gen_file, header = T) %>%
  mutate(u_ID = paste(aln, qlet, sep = "_")) %>%
  dplyr::select(u_ID, everything())

head(datafile_SNPs)

SNP_info <- read_csv("inputs/cad_gwas_snps.csv")
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

# Remove SNPs from SNP_info that weren't within the ALSPAC genotype data
SNP_info <- arrange(SNP_info, Lead_variant)
genotype_alleles <- arrange(genotype_alleles, SNP)
SNP_info[!(SNP_info$Lead_variant %in% genotype_alleles$SNP), "Lead_variant"]
SNP_info <- filter(SNP_info, Lead_variant %in% genotype_alleles$SNP)
stopifnot(SNP_info$Lead_variant == genotype_alleles$SNP)

# Extract SNPs where the effect allele does not equal the reference allele
rev_SNPs <- genotype_alleles %>%
  arrange(SNP) %>%
  mutate(Correct_ref = ifelse(ref_allele == SNP_info$A1, "T", "F")) %>%
  filter(Correct_ref == "F") %>%
  .[["SNP"]]

## NB: Weightings are based on the log(OR) of each variant with regards to CAD

datafile_SNPs_W <- datafile_SNPs

# Change the coding of effect alleles where needed and weight the SNPs
variants <- SNP_info$Lead_variant
for (i in variants) {
  # Changes allele codings
  if (i %in% rev_SNPs) {
    datafile_SNPs_W[[i]] <- 2 - datafile_SNPs_W[[i]]
  }
  # Weights the SNPs
  if (!(i %in% HMGCR_SNPs)) {
    weight <- subset(SNP_info, Lead_variant == i)
    new_col_name <- paste(i, "_w", sep = "")
    datafile_SNPs_W[[new_col_name]] <- datafile_SNPs_W[[i]] * weight[["logOR"]]
  }
}

# ------------------------------------------------------------------
# Merge the ages and and extract useful info
# ------------------------------------------------------------------
source("R/merge.R")
# merge 
d <- left_join(df_main, datafile_SNPs_W) %>%
  dplyr::select(u_ID, aln, qlet, age, everything()) %>%
  .[complete.cases(.), ]

d
nrow(d) # 5907 people

# For cohort characteristics
# temp <- dplyr::select(d, u_ID, aln, qlet, age)
# write.table(temp, file = "cad_individual_ids.txt", col.names = T, row.names = F, quote = F, sep = "\t")

# Extract SNP names  
SNPs <-  colnames(d)[str_detect(colnames(d), "[\\d]_w$")] 

# Extract lipoproteins
lipoproteins <- c(nr_mnames[grep("-V?[HIL]DL", nr_mnames)], nr_mnames[grep("IDL", nr_mnames)])
non_lipo <- nr_mnames[!(nr_mnames %in% lipoproteins)]

# Re-order the data by unique identifiers
d <- arrange(d, aln, qlet)
length(mnames) - length(nr_mnames) #81 ratios removed

# ------------------------------------------------------------------
# Generate biological subsets
# ------------------------------------------------------------------
# Generate the lipoprotein subsets 
Large_VLDL <- grep(paste(c("xxl-VLDL", "xl-VLDL", "l-VLDL", "^VLDL-"), collapse = "|"), nr_mnames, value = TRUE)
Atherogenic_non_LDL <- grep(paste(c("m-VLDL", "s-VLDL", "xs-VLDL", "IDL"), collapse = "|"), nr_mnames, value = TRUE)
LDL <- grep(paste(c("l-LDL", "m-LDL", "s-LDL"), collapse = "|"), nr_mnames, value = TRUE)
Very_large_HDL <- grep(c("xl-HDL"), nr_mnames, value = TRUE) 
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
Other <- c("Ace", "AcAce", "bOHBut", "Crea", "Alb", "Gp")
Other <- c(Other, diameter, cholesterol, glycerides, phospholipids, apolipo, fatty_acids, amino_acids,
           glycolysis)

# subsets <- list(Large_VLDL, Atherogenic_non_LDL, LDL, V_Large_HDL, Large_HDL, Small_HDL, diameter, cholesterol, glycerides, phospholipids, apolipo, fatty_acids, amino_acids, glycolysis, other)
subsets <- list(Large_VLDL, Atherogenic_non_LDL, LDL, Very_large_HDL, Large_HDL, Small_HDL, Other)
names(subsets) <- c("Large_VLDL", "Atherogenic_non_LDL", "LDL", "Very_large_HDL", "Large_HDL", "Small_HDL", "Other")

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

source("R/cluster_metabolites.R")

write.table(subset_df, file = "outputs/tables/subset_table.txt", quote = F, col.names = T, row.names = F, sep = "\t")

print("setup complete - please proceed to the CAD-GRS analysis")


