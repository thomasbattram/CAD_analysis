#  ----------------------------------------------------------------------------
# HMGCR SNP analysis 
# ----------------------------------------------------------------------------

HMGCR_SNP_results_nr <- lapply(HMGCR_SNPs, function(x) {linearRegress(x, nr_mnames, d, "age")})

names(HMGCR_SNP_results_nr) <- HMGCR_SNPs

# ----------------------------------------------------------------------------
# Extract the SNPs and metabolites with "significant" associations
# ----------------------------------------------------------------------------

HMGCR_SNPs_p_nr <- 0.05 / (length(HMGCR_SNPs) * length(nr_mnames))

# Add Bonferroni and FDR corrected p values to each of the results 
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
    } else if (type == 0.05) {
      res <- filter(data[[i]], `Pr(>|t|)` < 0.05)
    }
    if (nrow(res) > 0) {
      output[[i]] <- res
    }
  }
  output
}

# ----------------------------------------------------------------------------
# Make tables of results
# ----------------------------------------------------------------------------


sig_nr <- extract_sig_hits(HMGCR_SNP_results_nr)
for (i in names(sig_nr)) {
  sig_nr[[i]] <- sig_nr[[i]] %>%
    arrange(`Pr(>|t|)`) %>%
    mutate(P = make_pretty(`Pr(>|t|)`, 3)) %>%
    mutate(FDR = make_pretty(FDR, 3)) %>%
    mutate(low_CI = make_pretty(`2.5 %`, 3)) %>%
    mutate(up_CI = make_pretty(`97.5 %`, 3)) %>%
    mutate(`Estimate (95% CI)` = paste0(make_pretty(Estimate, 3), " (", low_CI, ", ", up_CI, ")")) %>%
    dplyr::select(Metabolite, `Estimate (95% CI)`, P, FDR)
}

sig_nr_FDR <- extract_sig_hits(HMGCR_SNP_results_nr, type = "fdr")
for (i in names(sig_nr_FDR)) {
  sig_nr_FDR[[i]] <- sig_nr_FDR[[i]] %>%
    arrange(`Pr(>|t|)`) %>%
    mutate(P = make_pretty(`Pr(>|t|)`, 3)) %>%
    mutate(FDR = make_pretty(FDR, 3)) %>%
    mutate(low_CI = make_pretty(`2.5 %`, 3)) %>%
    mutate(up_CI = make_pretty(`97.5 %`, 3)) %>%
    mutate(`Estimate (95% CI)` = paste0(make_pretty(Estimate, 3), " (", low_CI, ", ", up_CI, ")")) %>%
    dplyr::select(Metabolite, `Estimate (95% CI)`, P, FDR)
}

sig_nr_nom <- extract_sig_hits(HMGCR_SNP_results_nr, type = 0.05)
for (i in names(sig_nr_nom)) {
  sig_nr_nom[[i]] <- sig_nr_nom[[i]] %>%
    arrange(`Pr(>|t|)`) %>%
    mutate(P = make_pretty(`Pr(>|t|)`, 3)) %>%
    mutate(FDR = make_pretty(FDR, 3)) %>%
    mutate(low_CI = make_pretty(`2.5 %`, 3)) %>%
    mutate(up_CI = make_pretty(`97.5 %`, 3)) %>%
    mutate(`Estimate (95% CI)` = paste0(make_pretty(Estimate, 3), " (", low_CI, ", ", up_CI, ")")) %>%
    dplyr::select(Metabolite, `Estimate (95% CI)`, P, FDR)
}



write.table(sig_nr_nom[["rs17238484"]], file = "outputs/tables/rs17238484_sig_assoc.txt", quote = F, col.names = T, row.names = F, sep = "\t")
write.table(sig_nr_nom[["rs12916"]], file = "outputs/tables/rs12916_sig_assoc.txt", quote = F, col.names = T, row.names = F, sep = "\t")



