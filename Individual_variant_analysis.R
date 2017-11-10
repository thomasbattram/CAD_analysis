# ------------------------------------------------------------------
# Individual variant-metabolite analysis 
# ------------------------------------------------------------------

#drop the SNP with a log odds ratio of 0
SNPs <- SNPs %>%
  .[!(. %in% log_odd_SNP)]

d <- select(d, -one_of(log_odd_SNP))

indi_SNP_results_nr <- lapply(SNPs, function(x) {linearRegress(x, nr_mnames, d)})

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

write.table(names_sig_nr, paste0("outputs/other/", as.character(age), "/significant_SNPs.txt"), quote = F, col.names = F, row.names = F, sep = "\t")
write.table(names_sig_nr_FDR, paste0("outputs/other/", as.character(age), "/FDR_significant_SNPs.txt"), quote = F, col.names = F, row.names = F, sep = "\t")

workbook <- createWorkbook()
for (i in names_sig_nr_FDR) {
  temp_dat <- sig_nr_FDR[[i]] %>%
    select(Metabolite, everything()) %>%
    arrange(`Pr(>|t|)`)
  addWorksheet(wb = workbook, sheetName = i, gridLines = TRUE)
  writeDataTable(wb = workbook, sheet = i, x = temp_dat)
}
saveWorkbook(workbook, file = paste0("outputs/tables/", as.character(age), "/FDR_significant_SNP-metab_associations.xlsx"), overwrite = TRUE)

# ------------------------------------------------------------------
# Heatmap - beta-coef, data-driven clustering
# ------------------------------------------------------------------
source("R/Dendrogram_production.R")
load(file = "inputs/Pden_ColCol_variables_for_HeatMap.Rdata")


# Produce heatmaps using differing cluster methods and effect values
cluster_method <- c("data_driven", "biological")
data_type <- c("Pr(>|t|)", "Estimate")

for (i in 1:2) {
  data <- data_type[i]
  db <- sapply(indi_SNP_results_nr, function(x) {out = x[, data]; return(out)})
  rownames(db) <- rownames(indi_SNP_results_nr[[1]])
  db <- as.data.frame(db) %>%
    mutate(Metabolite = rownames(.)) %>%
    left_join(subset_df) %>%
    arrange(subset)

  rownames(db) <- db[["Metabolite"]]

  db <- select(db, -one_of(colnames(subset_df)))

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
      ColCol2[names(ColCol2) == "Other"] <- "white"
      Colv <- FALSE
      den <- "none"
      ColSC <- ColCol2
      fin_dat <- db
    }
    fin_dat <- as.matrix(fin_dat)

    pdf(paste0("outputs/heatmaps/", as.character(age), "/all_SNPs_vs_nr_metabs_", cluster, "_", data, "_heatmap.pdf"), width = 15, height = 10)
    heatmap.2( t(fin_dat), breaks = b, key = key, trace = "none", scale = "none", col = hmcol, rowsep = 1:62 , cexRow = 0.8, cexCol = 0.65, dendrogram = den , Colv =  Colv, Rowv = TRUE, ColSideColors = ColSC)
    dev.off()
  }
}


