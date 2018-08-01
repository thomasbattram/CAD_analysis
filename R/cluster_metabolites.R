# ------------------------------------------------------------------
# A first look at the NEO pilot data
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Read in the Data
# ------------------------------------------------------------------
dat <- df_main

# Remove individuals with missing data
dat <- dat %>%
	.[complete.cases(dat), ] %>%
	dplyr::select(u_ID, one_of(nr_mnames)) %>%
	.[, order(colnames(.))]
dim(dat)
colnames(dat)
# ------------------------------------------------------------------
#  Rank Normal Transformation of anthropetric 
#  and metabolite data
# ------------------------------------------------------------------

RNTdat <- apply(dplyr::select(dat, -u_ID), 2, rntransform)

# ------------------------------------------------------------------
#  Repeat analysis but use the Pearson Rho
#  values as distances for the dendrogram
# ------------------------------------------------------------------
D <-  as.dist(1 - abs(cor(RNTdat)))
PRhoTree <- hclust(D, method = "complete")
k <- cutree(PRhoTree, h = 0.2) ## distances are in 1 - Pearson Rho distances
table(k) ## 42 Clusters

# Plot dendrogram
pdf(paste0("outputs/other/Clustering_DendRect.pdf"), width = 20, height = 12)
par(mfrow = c(1,1))
plot(PRhoTree, hang = -1, cex = 0.5, main = "1 - Pearson Rho Cluster Dendrogram", col = "royalblue", lwd = 2)
rect.hclust(PRhoTree , h = 0.20, border = "red")
dev.off()

# NEW CODE USING GGPLOT!!!!
hcdata <- dendro_data(PRhoTree, type = "rectangle")
str(hcdata)
hc_clusters <- data.frame(label = names(k), clust = k)
rownames(hc_clusters) <- NULL

hc_labs <- label(hcdata) %>%
	left_join(hc_clusters) %>%
	left_join(subset_df, by = c("label" = "Metabolite")) %>%
	mutate(subset = ifelse(is.na(group), "Other", subset))

temp_hc_labs <- list()
# for (i in unique(hc_labs$clust)) {
# 	temp <- dplyr::filter(hc_labs, clust == i)
# 	temp2 <- list()
# 	for (j in unique(temp$subset)) {
# 		temp2 <- dplyr::filter(temp, subset == j)
# 		x <- temp2[["x"]]
# 		nr <- 1
# 		while (temp2[nr, "x"] - temp2[nr+1, "x"] == -1) {

# 		}
# 		temp2[[j]] <- dplyr::filter(temp, subset == j) %>%
# 			mutate(xmin = min(x) - 0.5) %>% # CHANGED FROM 0.4 TO 0.5!!!
# 			mutate(xmax = max(x) + 0.5)
# 	}
# 	temp_hc_labs[[i]] <- do.call(rbind, temp2)
# }

hc_labs %>%
	mutate(xmin = x - 0.5) %>%
	mutate(xmax = x + 0.5) -> hc_labs2


# hc_labs2 <- do.call(rbind, temp_hc_labs)
hc_labs2 <- arrange(hc_labs2, x)

p <- ggplot() + 
  geom_segment(data=segment(hcdata), aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_text(data=hc_labs2, aes(x=x, y=y, label=label, hjust=1), angle = 90, size=2) + 
  geom_hline(yintercept = 0.2, colour = "red") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
  		axis.title.x = element_blank(), legend.title = element_blank()) +
  geom_rect(data = hc_labs2, aes(xmin = xmin, xmax = xmax, ymin = -0.2, ymax = -0.065, fill = subset),
			size = 0.01) +
  labs(y = "Height", title = "1 - Pearson Rho Cluster Dendrogram")

p <- p + scale_fill_manual(breaks = c("Small_HDL", "Large_HDL", "Very_Large_HDL", "LDL", 
									  "Atherogenic_non_LDL", "Large_VLDL", "Other"), 
                           values=c("Small_HDL" = "red", "Large_HDL" = "blue", "Very_Large_HDL" = "green", "LDL" = "yellow", 
                           			"Atherogenic_non_LDL" = "purple", "Large_VLDL" = "orange", "Other" = "grey"))


ggsave("outputs/other/new_dendrogram.pdf", plot = p, width = 20, height = 12)

# ------------------------------------------------------------------
# Produce dendrogram and colour variables for heatmaps
# ------------------------------------------------------------------

Pden <- as.dendrogram(PRhoTree)

M <- RNTdat
PMat <- as.dist(1 - abs(cor(M, use = "p")))
PMat_Tree <- hclust(PMat, method = "complete")
K <- cutree(PMat_Tree, h = 0.2)
length(unique(K))
temp <- which(table(K) == 1); K[K %in% temp] = 0; length(unique(K))
x <- as.numeric(names(table(K))) + 100
K <- K + 100
for(i in 1:length(x)){
  r  <- which(K == x[i])
  K[r] <- i
}

CC <- sample(colorRampPalette(brewer.pal(9, "Set1"))(length(unique(K))))
names(CC) = 1:length(unique(K))
ColCol = CC[K]

save(Pden, ColCol, file = "inputs/Pden_ColCol_variables_for_HeatMap.Rdata")

# ------------------------------------------------------------------
# Produce a colour variable based on biological grouping
# ------------------------------------------------------------------

CC2 <- sample(colorRampPalette(brewer.pal(9, "Set1"))(length(subsets)))
names(CC2) <- c(names(subsets))

ordered_df <- arrange(subset_df, subset)

K2 <- ordered_df[["Metabolite"]]
names(K2) <- ordered_df[["subset"]]

ColCol2 <- CC2[names(K2)]

# ------------------------------------------------------------------
# Produce a new colour variable 
# ------------------------------------------------------------------

# colours = biological groups
# order = alphabetical metabolites

ordered_df2 <- arrange(subset_df, Metabolite)
K3 <- ordered_df2[["Metabolite"]]
names(K3) <- ordered_df2[["subset"]]

ColCol3 <- CC2[names(K3)]

## add clusters to subset_df
clusters <- data.frame(Metabolite = names(k), cluster = k)
rownames(clusters) <- NULL

subset_df <- subset_df %>%
	left_join(clusters)

dplyr::filter(subset_df, !is.na(group))
