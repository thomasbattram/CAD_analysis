##############################################
##
##  A first look at the NEO pilot data
##
##############################################

###############################
## Read in the Data
###############################
dat <- df_main

## Remove individuals with missing data
dat <- dat %>%
	.[complete.cases(dat), ] %>%
	dplyr::select(u_ID, one_of(nr_mnames)) %>%
	.[, order(colnames(.))]
dim(dat)
colnames(dat)
##########################################
##  Rank Normal Transformation of anthropetric 
##  and metabolite data
##########################################

#rank normal transformation - produces perfectly normal distribution of some data
RNTdat <- apply(dplyr::select(dat, -u_ID), 2, rntransform)

############################################
##
##  Repeat analysis but use the Pearson Rho
##  values as distances for the dendrogram
##
############################################
#M = RNTdat
D <-  as.dist( 1 - abs(cor(RNTdat)) )
PRhoTree <- hclust(D, method = "complete")
k <- cutree(PRhoTree, h = 0.2) ## distances are in 1 - Pearson Rho distances CHANGE HERE TO GET A PLOT THAT FITS THE HEATMAPS BETTER
table(k) ## 43 Clusters

## Plot dendrogram
pdf("outputs/Clustering_DendRect.pdf", width = 20, height = 12)
par(mfrow = c(1,1))
plot(PRhoTree, hang = -1, cex = 0.5, main = "1 - Pearson Rho Cluster Dendrogram", col = "royalblue", lwd = 2)
rect.hclust( PRhoTree , h = 0.20, border = "red")
dev.off()

############################

Pden <- as.dendrogram(PRhoTree)

M <- RNTdat
PMat <- as.dist(1 - abs(cor(M, use="p")))
PMat_Tree <- hclust(PMat , method="complete")
K <- cutree(PMat_Tree, h = 0.2)
length(unique(K))
temp <- which(table(K) == 1); K[ K %in% temp ] = 0; length(unique(K))
x <- as.numeric(names(table(K))) + 100
K <- K + 100
for(i in 1:length(x)){
  r  <- which(K == x[i])
  K[r] <- i
}
##
CC <- sample( colorRampPalette(brewer.pal(9,"Set1"))(length(unique(K))) )
names(CC) = 1:length(unique(K))
ColCol = CC[K]
##
save(Pden, ColCol, file = "inputs/Pden_ColCol_variables_for_HeatMap.Rdata")


################# ---------------------------------------------------------------------
#
# 	Experimenting with grouping based on lipoproteins
#
################# ---------------------------------------------------------------------

subsets
non_lipo <- nr_mnames[!(nr_mnames %in% lipoproteins)]

#new sets
diameter <- non_lipo[grep("-D", non_lipo)]
cholesterol <- non_lipo[grep("-C$", non_lipo)]
glycerides <- non_lipo[grep("-TG$", non_lipo)]
phospholipids <- c("DAG", "Tot-PG", "PC", "Tot-Cho")
apolipo <- non_lipo[grep("Apo", non_lipo)]
fatty_acids <- c("Tot-FA", "FALen", "UnsatDeg", "DHA", "LA", "CLA", "FAw3", "FAw6", "PUFA", "MUFA", "SFA")
amino_acids <- c("Ala", "Gln", "His", "Ile", "Leu", "Val", "Phe", "Tyr")
glycolysis <- c("Glc", "Lac", "Pyr", "Cit", "glucose", "insulin")
other <- c("Ace", "AcAce", "bOHBut", "Crea", "Alb", "Gp")



new_sets <- list(diameter, cholesterol, glycerides, phospholipids, apolipo, fatty_acids, amino_acids, glycolysis, other)
names(new_sets) <- c("diameter", "cholesterol", "glycerides", "phospholipids", "apolipo", "fatty_acids", "amino_acids", "glycolysis", "other")

test <- unlist(new_sets)

names(test) <- str_replace(names(test), "[0-9]", "")

CC2 <- sample( colorRampPalette(brewer.pal(9,"Set1"))(length(subsets) + length(new_sets) ))
names(CC2) <- c(names(subsets), names(new_sets))
###DOESN'T WORK - CHANGE IT HERE!!!!!

temp_dat <- CAD_score_lr_nr %>%
	mutate(new_subset = ifelse(is.na(subset), case_when(
		Metabolite %in% diameter ~ "diameter",
		Metabolite %in% cholesterol ~ "cholesterol",
		Metabolite %in% glycerides ~ "glycerides",
		Metabolite %in% phospholipids ~ "phospholipids",
		Metabolite %in% apolipo ~ "apolipoprotein",
		Metabolite %in% fatty_acids ~ "fatty_acids",
		Metabolite %in% amino_acids ~ "amino_acids",
		Metabolite %in% glycolysis ~ "glycolysis",
		Metabolite %in% other ~ "other"), as.character(subset))) %>%
	arrange(new_subset)
K2 <- temp_dat[["Metabolite"]]
names(K2) <- temp_dat[["new_subset"]]

ColCol2 <- CC2[names(K2)]





