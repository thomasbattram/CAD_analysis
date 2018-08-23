# ------------------------------------------------------------------
# Individual variant-metabolite analysis 
# ------------------------------------------------------------------

indi_SNP_results_nr <- lapply(SNPs, function(x) {linearRegress(x, nr_mnames, d, "age")})
HMGCR_SNP_results_nr <- lapply(HMGCR_SNPs, function(x) {linearRegress(x, nr_mnames, d, "age")})
names(HMGCR_SNP_results_nr) <- HMGCR_SNPs

SNP_names <- unlist(strsplit(SNPs, "_w"))

names(indi_SNP_results_nr) <- SNP_names

# ------------------------------------------------------------------
# Extract the SNPs and metabolites with "significant" associations 
# ------------------------------------------------------------------

indi_SNPs_p_nr <- 0.05 / (length(SNPs) * length(nr_mnames))

# Add Bonferroni and FDR corrected p values to each of the results 
for(i in 1:length(indi_SNP_results_nr)) {
  indi_SNP_results_nr[[i]]$Bonferroni <- p.adjust(indi_SNP_results_nr[[i]]$`Pr(>|t|)`, method = "bonferroni", n = (length(SNPs) * length(nr_mnames)))
  indi_SNP_results_nr[[i]]$FDR <- p.adjust(indi_SNP_results_nr[[i]]$`Pr(>|t|)`, method = "fdr", n = (length(SNPs) * length(nr_mnames)))
}

save(indi_SNP_results_nr, file = "outputs/other/all_indi_res.RData")

p_score_nr <- 0.05/length(nr_mnames)

# Make a list of the results containing the "significant" p values
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

sig_nr <- extract_sig_hits(indi_SNP_results_nr)
sig_nr_FDR <- extract_sig_hits(indi_SNP_results_nr, type = "fdr")

# Write the results into tables
names_sig_nr <- names(sig_nr)
names_sig_nr_FDR <- names(sig_nr_FDR)

write.table(names_sig_nr, "outputs/other/significant_SNPs.txt", quote = F, col.names = F, row.names = F, sep = "\t")
write.table(names_sig_nr_FDR, "outputs/other/FDR_significant_SNPs.txt", quote = F, col.names = F, row.names = F, sep = "\t")

CAD_tab <- arrange(CAD_score_lr_nr, `Pr(>|t|)`) %>%
  mutate(P = make_pretty(`Pr(>|t|)`, 3)) %>%
  mutate(low_CI = make_pretty(`2.5 %`, 3)) %>%
  mutate(up_CI = make_pretty(`97.5 %`, 3)) %>%
  mutate(`Estimate (95% CI)` = paste0(make_pretty(Estimate, 3), " (", low_CI, ", ", up_CI, ")")) %>%
  dplyr::select(Metabolite, `Estimate (95% CI)`, P)

workbook <- createWorkbook()
for (i in names_sig_nr_FDR) {
  temp_dat <- sig_nr_FDR[[i]] %>%
    mutate(P = make_pretty(`Pr(>|t|)`, 3)) %>%
    mutate(FDR = make_pretty(FDR, 3)) %>%
    mutate(low_CI = make_pretty(`2.5 %`, 3)) %>%
    mutate(up_CI = make_pretty(`97.5 %`, 3)) %>%
    mutate(`Estimate (95% CI)` = paste0(make_pretty(Estimate, 3), " (", low_CI, ", ", up_CI, ")")) %>%
    dplyr::select(Metabolite, `Estimate (95% CI)`, P, FDR) %>%
    arrange(P)
  addWorksheet(wb = workbook, sheetName = i, gridLines = TRUE)
  writeDataTable(wb = workbook, sheet = i, x = temp_dat)
}
saveWorkbook(workbook, file = "outputs/tables/FDR_significant_SNP-metab_associations.xlsx", overwrite = TRUE)

# ------------------------------------------------------------------
# Heatmaps
# ------------------------------------------------------------------
load(file = "inputs/Pden_ColCol_variables_for_HeatMap.Rdata")

# Produce heatmaps using differing cluster methods and effect values
cluster_method <- c("data_driven", "biological")
data_type <- c("Pr(>|t|)", "Estimate")
i=1
j=2
heat_dat <- c(indi_SNP_results_nr, HMGCR_SNP_results_nr)
for (i in 1:2) {
  data <- data_type[i]
  # Extract the data needed and arrange the metabolites
  db <- sapply(heat_dat, function(x) {out = x[, data]; return(out)})
  rownames(db) <- rownames(heat_dat[[1]])
  db <- as.data.frame(db) %>%
    mutate(Metabolite = rownames(.)) %>%
    left_join(subset_df) %>%
    arrange(subset)

  rownames(db) <- db[["Metabolite"]]

  db <- dplyr::select(db, -one_of(colnames(subset_df)))


  # Extract variables required for heatmap function
  if (data == "Estimate") {
    b <- seq(from = -5, to = 5, by = 1)
    hmcol <- colorRampPalette(brewer.pal(9, "RdBu"))
    key <- TRUE
  } else if (data == "Pr(>|t|)") {
    b <- c(0, 1e-8, 1e-6, 1e-4, 1e-2, 0.05, .25, 0.5, 0.75, 1) 
    hmcol <- brewer.pal(11, "Spectral")[1:5]
    sig <- brewer.pal(11, "Spectral")[1:5]
    hmcol <- c(sig,"grey90","grey75","grey50", "grey45")
    key <- FALSE
  }
  for (j in 1:2) {
    cluster <- cluster_method[j]
    if (cluster == "data_driven") {
      ColCol[names(ColCol) %in% 1] <- "white"
      Colv <- Pden
      den <- "both"
      ColSC <- ColCol
      fin_dat <- db[order(row.names(db)), ]
    } else if (cluster == "biological") {
      Colv <- FALSE
      den <- "none"
      ColSC <- ColCol2
      fin_dat <- db
    }

    fin_dat <- as.matrix(fin_dat)

    write.table(fin_dat, file = paste0("outputs/heatmaps/", cluster, "_", data, "dat.txt"), row.names = T, col.names = T, quote = F, sep = "\t")

    # Remove HMGCR SNPs from the original heatmap
    heat1 <- fin_dat[, !(colnames(fin_dat) %in% HMGCR_SNPs)]

    pdf(paste0("outputs/heatmaps/all_SNPs_vs_nr_metabs_", cluster, "_", data, "_heatmap.pdf"), width = 15, height = 10)
    heatmap <- heatmap.2( t(heat1), breaks = b, key = key, trace = "none", scale = "none", col = hmcol, rowsep = 1:ncol(fin_dat) , cexRow = 0.5, cexCol = 0.65, dendrogram = "none", Colv =  Colv, Rowv = TRUE, ColSideColors = ColSC, margins =c(5,9))
    print(heatmap)
    dev.off()
    ####### The colours aren't correct - need to find a way to change the ordering step!!!!
    # Make a heatmap with just the SNPs associated with metabolites and the HMGCR SNPs
    new_fin_dat <- fin_dat[, c(names_sig_nr_FDR, HMGCR_SNPs)]

    gene_info <- SNP_info %>%
      dplyr::select(Lead_variant, gene) %>%
      filter(Lead_variant %in% colnames(new_fin_dat))

    index <- match(colnames(new_fin_dat), gene_info$Lead_variant)
    gene_info <- gene_info[index, ]

    gene_info[gene_info$gene == "BUD13_ZNF259_APOA5", "gene"] <- "APOA5"
    gene_info[gene_info$gene == "TOMM40_APOE_APOC1", "gene"] <- "APOE_APOC1"

    colnames(new_fin_dat) <- paste0(colnames(new_fin_dat), "\n", gene_info$gene)


    pdf(paste0("outputs/heatmaps/sig_SNPs_vs_nr_metabs_", cluster, "_", data, "_heatmap.pdf"), width = 15, height = 10)
    heatmap.2( t(new_fin_dat), breaks = b, key = key, trace = "none", scale = "none", col = hmcol, rowsep = 1:ncol(new_fin_dat) , cexRow = 1.15, cexCol = 0.65, dendrogram = den , Colv =  Colv, Rowv = TRUE, ColSideColors = ColSC, margins =c(5,9))
    dev.off()

  }
}

