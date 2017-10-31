##########################################
# Individual variant-metabolite analysis #
##########################################

#source the linear regression function - mostly made by Qin Wang
source("R/Linear_regression_func.R")

#drop the SNP with a log odds ratio of 0
SNPs <- SNPs %>%
  .[!(. %in% log_odd_SNP)]

d <- select(d, -one_of(log_odd_SNP))

#
indi_SNP_results_nr <- lapply(SNPs, function(x) {linearRegress(x, nr_mnames, d)})

SNP_names <- unlist(strsplit(SNPs, "_w"))

names(indi_SNP_results_nr) <- SNP_names

############################################################################
###6. Extract the SNPs and metabolites with significant associations #######
############################################################################

indi_SNPs_p_nr <- 0.05 / (length(SNPs) * length(nr_mnames))

#Add Bonferroni and FDR corrected p values to each of the results 
for(i in 1:length(indi_SNP_results_nr)) {
  indi_SNP_results_nr[[i]]$Bonferroni <- p.adjust(indi_SNP_results_nr[[i]]$`Pr(>|t|)`, method = "bonferroni", n = (length(SNPs) * length(nr_mnames)))
  indi_SNP_results_nr[[i]]$FDR <- p.adjust(indi_SNP_results_nr[[i]]$`Pr(>|t|)`, method = "fdr", n = (length(SNPs) * length(nr_mnames)))
}

p_score_nr <- 0.05/length(nr_mnames)

#make a list of the results containing the significant p values
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

#############################################################################
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

#######################################
#### HEATMAPS - beta-coefs #######
#######################################
###If using Pden and ColCol:
source("R/Dendrogram_production.R")
############

###Extract the beta-coefficients 

db <- sapply(indi_SNP_results_nr, function(x) {out = x[, 1]; return(out)})

rownames(db) <- rownames(indi_SNP_results_nr[[1]])


#put the metabolites in order
db <- as.data.frame(db)


db <- db[order(row.names(db)), ]

db <- as.matrix(db)

#sum(rownames(db) == K2)
###Plot the heatmaps

#Load Pre-computed Plotting variables
#(Metabolite Clusters From ALSPAC Data and Metabolite Cluster ColorScheme)
load( file = "inputs/Pden_ColCol_variables_for_HeatMap.Rdata")
ColCol[names(ColCol) %in% 1] <- "white"

#
ColCol2[names(ColCol2) == "Other"] <- "white"
#

#for the colour key scale - allows clearer visualisation of associations
b <- seq(from = -5, to = 5, by = 1)

### Plotting Colors
hmcol <- colorRampPalette(brewer.pal(9, "RdBu"))


#all SNPs vs no ratios metabs 
pdf(paste0("outputs/heatmaps/", as.character(age), "/all_SNPs_vs_nr_metabs_coef.pdf"), width = 15, height = 10)
heatmap.2( t(db), breaks = b, trace = "none", scale = "none", col = hmcol, rowsep = 1:62 , cexRow = 0.8, cexCol = 0.65, dendrogram = "both" , Colv =  Pden, Rowv = TRUE, ColSideColors = ColCol)
dev.off()


###########################
###	HEATMAPS - P values ###
###########################

db <- sapply(indi_SNP_results_nr, function(x) {out = x[, 4]; return(out)})

rownames(db) <- rownames(indi_SNP_results_nr[[1]])


SNP_loci <- SNP_info %>%
  filter(Lead_variant %in% colnames(db))

SNP_loci <- paste(SNP_loci$Lead_variant, SNP_loci$Locus, sep = "_")



#put the metabolites in order
db <- as.data.frame(db)

db <- select(db, sort(colnames(db)))
colnames(db) <- SNP_loci

db <- db[order(row.names(db)), ]

db <- as.matrix(db)
head(db)
#Set groups into different variables

###Plot the heatmaps

#Load Pre-computed Plotting variables
#(Metabolite Clusters From ALSPAC Data and Metabolite Cluster ColorScheme)
load( file = "inputs/Pden_ColCol_variables_for_HeatMap.Rdata")
ColCol[names(ColCol) %in% 1] <- "white"
 
b <- c(0, 1e-8, 1e-6, 1e-4, 1e-2, 0.05, .25, 0.5, 0.75, 1)

### Plotting Colors
hmcol <- brewer.pal(11, "Spectral")[1:5]
sig <- brewer.pal(11, "Spectral")[1:5]
hmcol <- c(sig,"grey90","grey75","grey50", "grey45")

## p value PLOT
#all SNPs vs no ratios metabs
pdf(paste0("outputs/heatmaps/", as.character(age), "/all_SNPs_vs_nr_metabs_pval.pdf"), width = 15, height = 10)
heatmap.2( t(db), key = FALSE, breaks = b,  trace = "none", scale = "none", col = hmcol, rowsep = 1:62 , cexRow = 0.8, cexCol = 0.65, dendrogram = "column" , Colv =  Pden, Rowv = TRUE, ColSideColors = ColCol)
dev.off()



print("Individual variant analysis complete - please proceed to the New-GRS analysis")

########### ----------------------------------------------------------------------------
#
# Experimenting with changing heatmap order
#
########### ----------------------------------------------------------------------------

#Coefs

db <- sapply(indi_SNP_results_nr, function(x) {out = x[, 1]; return(out)})

rownames(db) <- rownames(indi_SNP_results_nr[[1]])

#
stopifnot(sum(temp_dat$Metabolite %in% rownames(db)) == length(nr_mnames))
#

#
db <- as.data.frame(db) %>%
  mutate(Metabolite = rownames(.)) %>%
  left_join(temp_dat) %>%
  arrange(new_subset)

rownames(db) <- db[["Metabolite"]]

db <- select(db, -one_of(colnames(temp_dat)))

db <- as.matrix(db)
#

#for the colour key scale - allows clearer visualisation of associations
b <- seq(from = -5, to = 5, by = 1)

### Plotting Colors
hmcol <- colorRampPalette(brewer.pal(9, "RdBu"))

#
pdf(paste0("outputs/heatmaps/", as.character(age), "/all_SNPs_vs_nr_metabs_coef_TEST.pdf"), width = 15, height = 10)
heatmap.2( t(db), breaks = b, trace = "none", scale = "none", col = hmcol, rowsep = 1:62 , cexRow = 0.8, cexCol = 0.65, dendrogram = "row", Colv = FALSE, Rowv = TRUE, ColSideColors = ColCol2)
dev.off()
#


#P vals

db <- sapply(indi_SNP_results_nr, function(x) {out = x[, 4]; return(out)})

rownames(db) <- rownames(indi_SNP_results_nr[[1]])



#
stopifnot(sum(temp_dat$Metabolite %in% rownames(db)) == length(nr_mnames))
#

#
db <- as.data.frame(db) %>%
  select(sort(colnames(db)))

#Sanity check
stopifnot(colnames(db) == unlist(str_split(SNP_loci, "_"))[seq(1, 110, by = 2)])

#Add in gene names
colnames(db) <- SNP_loci


db <- mutate(db, Metabolite = rownames(db)) %>%
  left_join(temp_dat) %>%
  arrange(new_subset)
colnames(db)


rownames(db) <- db[["Metabolite"]]
colnames(db)
db <- select(db, -one_of(colnames(temp_dat)))
colnames(db)
db <- as.matrix(db)
#

#breaks
b <- c(0, 1e-8, 1e-6, 1e-4, 1e-2, 0.05, .25, 0.5, 0.75, 1)

### Plotting Colors
hmcol <- brewer.pal(11, "Spectral")[1:5]
sig <- brewer.pal(11, "Spectral")[1:5]
hmcol <- c(sig,"grey90","grey75","grey50", "grey45")

#
pdf(paste0("outputs/heatmaps/", as.character(age), "/all_SNPs_vs_nr_metabs_pval_TEST.pdf"), width = 15, height = 10)
heatmap.2( t(db), key = FALSE, breaks = b, trace = "none", scale = "none", col = hmcol, rowsep = 1:62 , cexRow = 0.8, cexCol = 0.65, dendrogram = "none", Colv = FALSE, Rowv = TRUE, ColSideColors = ColCol2, lwid = c(0.1, 4), lhei = c(0.1, 4, 1))
dev.off()
#









