######################
# HMGCR SNP analysis #
######################


#source the linear regression function - mostly made by Qin Wang
source("R/Linear_regression_func.R")

HMGCR_SNP_results_nr <- lapply(HMGCR_SNPs, function(x) {linearRegress(x, nr_mnames, d)})

names(HMGCR_SNP_results_nr) <- HMGCR_SNPs

############################################################################
### Extract the SNPs and metabolites with significant associations #######
############################################################################

HMGCR_SNPs_p_nr <- 0.05 / (length(HMGCR_SNPs) * length(nr_mnames))

#Add Bonferroni and FDR corrected p values to each of the results 
for(i in 1:length(HMGCR_SNP_results_nr)) {
  HMGCR_SNP_results_nr[[i]]$Bonferroni <- p.adjust(HMGCR_SNP_results_nr[[i]]$`Pr(>|t|)`, method = "bonferroni", n = (length(HMGCR_SNPs) * length(nr_mnames)))
  HMGCR_SNP_results_nr[[i]]$FDR <- p.adjust(HMGCR_SNP_results_nr[[i]]$`Pr(>|t|)`, method = "fdr", n = (length(HMGCR_SNPs) * length(nr_mnames)))
}


extract_sig_hits <- function(data, type = "bon") {
  output <- list()
  for (i in names(data)) {
    if (type == "bon") {
      res <- filter(data[[i]], Bonferroni < 0.05)
    } else if (type == "fdr") {
      res <- filter(data[[i]], FDR < 0.05)
    }
    if (nrow(res) > 0) {
      output[[i]] <- res
    }
  }
  output
}

sig_nr <- extract_sig_hits(HMGCR_SNP_results_nr)
sig_nr_FDR <- extract_sig_hits(HMGCR_SNP_results_nr, type = "fdr")

##########################
# Make tables of results #
##########################

write.table(sig_nr_FDR[["rs17238484"]], file = paste0("outputs/other/", as.character(age), "/rs17238484_sig_assoc.txt"), quote = F, col.names = T, row.names = F, sep = "\t")
write.table(sig_nr_FDR[["rs12916"]], file = paste0("outputs/other/", as.character(age), "/rs12916_sig_assoc.txt"), quote = F, col.names = T, row.names = F, sep = "\t")








