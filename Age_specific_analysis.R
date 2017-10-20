###
#Age specific stuff
###
rm(list = ls())
library(cluster)
library(RColorBrewer)
library(gplots)
library(ape)
library(readstata13)
library(foreign)
library(xlsx)
library(GenABEL)
library(gtools)
library(ggplot2)
library(tidyverse)
library(haven)
library(stringr)
library(zeallot)

setwd("CAD_adolescent_analysis")

#the original data file with metabs
datafile_metabs <- read_dta("inputs/metabolite_data.dta") #was export_metabolomics

#Add in a unique identifier
datafile_metabs <- datafile_metabs %>%
  mutate(u_ID = paste(cidB9999, qlet, sep = "_")) %>%
  dplyr::select(u_ID, everything())

#### Need to alter metabolite names so they are paper presentable

#Name metabolites that are named OK in the metabolite dataset
fine_mets <- c("DAG", "PC", "ApoA1", "ApoB", "FALen", "UnsatDeg", "DHA", "LA", "CLA", "FAw3", "FAw6", "PUFA", "MUFA", "SFA", "Glc", "Lac", "Pyr", "Cit", "Ala", "Gln", "His", "Ile", "Leu", "Val", "Phe", "Tyr", "Ace", "AcAce", "bOHBut", "Crea", "Alb", "Gp", "glucose", "insulin")
#Name the lipoproteins that don't have a size suffix
diff_lipos <- c("VLDLD", "LDLD", "HDLD", "VLDLC", "LDLC", "HDLC", "HDL2C", "HDL3C", "VLDLTG", "LDLTG", "HDLTG")
new_diff_lipos <- paste(
  str_extract(diff_lipos, "V?[HIL]DL"),
  str_sub(str_extract(diff_lipos, "DL.+"), 3, 10),
  sep = "-"
  )
#Put together a vector of how the wrongly named non-lipo-protein metabolites should be named - in order of the metabolites in the column names of the df
new_non_lipos <- c("Serum-C", "Remnant-C", "Est-C", "Free-C", "Serum-TG", "DAG/TG", "Tot-PG", "TG/PG", "Tot-Cho", "ApoB/ApoA1", "Tot-FA", "DHA/FA", "LA/FA", "CLA/FA", "FAw3/FA", "FAw6/FA", "PUFA/FA", "MUFA/FA", "SFA/FA")
df_list <- list()
#F7 = age 7
#TF3 = age 15
#TF4 = age 17
suffixes <- c("_F7", "_TF3", "_TF4")
#i = suffixes[2]
#Make a loop to extract data at all 3 ages and change the metabolite names
for (i in suffixes) {
  df <- dplyr::select(datafile_metabs, u_ID, matches(paste(i, "$", sep = ""), ignore.case = FALSE))
  colnames(df) <- str_replace(colnames(df), i, "")
  lipos <- colnames(df)[matches("[HIL]DL", vars = colnames(df))] %>%
    .[!(. %in% diff_lipos)]
  sizes <- str_c(c("X+L", "L", "M" ,"S", "X+S"), collapse = "|")
  new_lipos <- paste(
    str_to_lower(str_extract(lipos, sizes)), #First part
    str_extract(lipos, "V?[HIL]DL"), #Second part
    str_sub(str_extract(lipos, "DL.+"), 3, 10), #Third part
    sep = "-"
    ) 
  colnames(df)[colnames(df) %in% diff_lipos] <- new_diff_lipos
  colnames(df)[colnames(df) %in% lipos] <- new_lipos
  non_lipos <- colnames(df)[!(colnames(df) %in% new_lipos | colnames(df) %in% new_diff_lipos)] %>%
    .[-1] %>%
    .[!(. %in% fine_mets)]

  colnames(df)[colnames(df) %in% non_lipos] <- new_non_lipos

  df_list[[i]] <- df
}

#Make 3 data frames that have only the complete cases for each age
df_f7 <- df_list[[1]] %>%
  .[complete.cases(.), ]
df_f7[["age"]] <- 7


df_tf3 <- df_list[[2]] %>%
  .[complete.cases(.), ] %>%
  select(-glucose, -insulin) #Remove glucose and insulin to increase complete cases sample size
df_tf3[["age"]] <- 15

df_tf4 <- df_list[[3]] %>%
  .[complete.cases(.), ] 
df_tf4[["age"]] <- 17

#Bind the data frames together
df_all <- rbind(df_f7, df_tf3, df_tf4)




#####Tukey function

Tukey <- function(dataframe, col_name) {
  three_IQR <- 3*IQR(dataframe[[col_name]], na.rm=T)
  quartiles <- summary(dataframe[[col_name]])[c(2,5)]
  lo_lim <- as.numeric(quartiles[1] - three_IQR)
  up_lim <- as.numeric(quartiles[2] + three_IQR)

  large_vals <- dataframe[dataframe[[col_name]] > lo_lim, ]
  vals <- large_vals[large_vals[[col_name]] < up_lim, ]

  adj_sites <- which(rownames(dataframe) %in% rownames(vals))
  dataframe <- dataframe[adj_sites, ]
}

#Trimming metabolites
tukey_res <- df_all
for (i in mnames) {
  tukey_res <- Tukey(tukey_res, i)
}
nrow(tukey_res) #only 2424 observations left!!!


mnames <- colnames(df_all)[-c(1, 232)]
nr_mnames <- mnames[-grep("_P|/", mnames)]

length(unique(df_all[["u_ID"]])) #Number of unique IDs (people) = 6990
nrow(df_all) #Total number of IDs = 11344
u_ID <- unique(df_all[["u_ID"]]) #Unique IDs
tab_ID <- table(df_all[["u_ID"]])
nrow(df_f7) #5353
nrow(df_tf3) #2937
nrow(df_tf4) #3054

#Change the data frames so that there is no cross over of IDs between the ages
df_f7_filt <- df_f7 %>%
  filter(!(u_ID %in% df_tf3[["u_ID"]])) %>%
  filter(!(u_ID %in% df_tf4[["u_ID"]]))

df_tf3_filt <- df_tf3 %>%
  filter(!(u_ID %in% df_f7[["u_ID"]])) %>%
  filter(!(u_ID %in% df_tf4[["u_ID"]]))

df_tf4_filt <- df_tf4 %>%
  filter(!(u_ID %in% df_f7[["u_ID"]])) %>%
  filter(!(u_ID %in% df_tf3[["u_ID"]]))

length(unique(df_all[["u_ID"]])) - nrow(df_f7_filt) - nrow(df_tf4_filt) - nrow(df_tf3_filt) #Number of IDs in all 3 datasets = 3162

#Function takes 2 datasets and produces a single dataset containing the same number of observations from both the input datasets
sample_random_IDs <- function(df_large, df_small, col_ID, col_diff = FALSE) {
  if (nrow(df_large) < nrow(df_small)) stop("df_large is smaller than df_small")
  if (col_diff != FALSE) {
    df_l_length <- nrow(df_large) / length(unique(df_large[[col_diff]]))
    df_s_length <- nrow(df_small) / length(unique(df_small[[col_diff]]))
    length_diff <-  df_l_length - df_s_length
    n_diff <- length(unique(df_large[[col_diff]])) + length(unique(df_small[[col_diff]]))
    } else {
      length_diff <- nrow(df_large) - nrow(df_small)
      n_diff <- 2
    }

  same_df <- df_large[df_large[[col_ID]] %in% df_small[[col_ID]], ] %>%
    sample_n((nrow(.) + length_diff) / n_diff)
  df1 <- df_large[!(df_large[[col_ID]] %in% same_df[[col_ID]]), ]
  df2 <- df_small[!(df_small[[col_ID]] %in% df1[[col_ID]]), ]
  if (length_diff %% 2 != 0) {
    df1 <- sample_n(df1, nrow(df1) - 1)
  }
  final_df <- rbind(df1, df2)
  final_df
}


#Set the seed so the results are reproducible
set.seed(2)
df_temp <- sample_random_IDs(df_tf4, df_tf3, 1)

#Check there are the same number of individuals at each age
table(df_temp[["age"]])
nrow(df_temp) / 2 #Number of each age group = 2126

set.seed(22)
df_f7_samp <- df_f7_filt %>%
  sample_n(nrow(df_temp) / 2) 

df_fin <- rbind(df_temp, df_f7_samp)

table(df_fin[["age"]])

###############################################
# Testing differences using statistical tests #
###############################################
df_fin_7 <- filter(df_fin, age == 7)
df_fin_17 <- filter(df_fin, age == 17)
df_fin_15 <- filter(df_fin, age == 15)
#Wilcoxon ranked sum test 1 - testing diff between ages 7 and 17
w_test_p <- c()
for (i in mnames) {
  x <- wilcox.test(df_fin_7[[i]], df_fin_17[[i]], alternative = "two.sided")
  p_val <- x[[3]]
  names(p_val) <- i
  w_test_p <- c(w_test_p, p_val)
}

k_test_p <- c()
for(i in mnames) {
  fom <- formula(paste("`", i, "`", " ~ ", "age", sep = ""))
  x <- kruskal.test(fom, data = df_fin)
  p_val <- x[[3]]
  names(p_val) <- i
  k_test_p <- c(k_test_p, p_val)
}

length(k_test_p)
length(k_test_p[k_test_p < 0.05]) #224 of 230 < 0.05

#For KW test:
# H0 = medians of all groups are equal
# H1 = median of at least one group is different from the median of one other group


#Wilcoxon ranked sum test 2 - testing diff between ages 15 and 17
w_test_p_2 <- c()
for (i in mnames) {
  x <- wilcox.test(df_fin_15[[i]], df_fin_17[[i]], alternative = "two.sided")
  p_val <- x[[3]]
  names(p_val) <- i
  w_test_p_2 <- c(w_test_p_2, p_val)
}

length(w_test_p) #230
length(w_test_p[w_test_p < 0.05]) #215 of 230 < 0.05
length(w_test_p_2[w_test_p_2 < 0.05]) #168 of 230 < 0.05

####### Using rank normalised data and t.test - p vals are all ~1 because of the rank normalisation...
df_fin_7_rnt <- df_fin_7
df_fin_17_rnt <- df_fin_17
df_fin_15_rnt <- df_fin_15
t_test_p <- c()
for (i in colnames(df_fin)[-c(1,232)]) {
  df_fin_7_rnt[[i]] <- as.numeric(df_fin_7_rnt[[i]])
  df_fin_7_rnt[[i]] <- rntransform(formula = df_fin_7_rnt[[i]], data = df_fin_7_rnt, family = gaussian)

  df_fin_17_rnt[[i]] <- as.numeric(df_fin_17_rnt[[i]])
  df_fin_17_rnt[[i]] <- rntransform(formula = df_fin_17_rnt[[i]], data = df_fin_17_rnt, family = gaussian)

  df_fin_15_rnt[[i]] <- as.numeric(df_fin_15_rnt[[i]])
  df_fin_15_rnt[[i]] <- rntransform(formula = df_fin_15_rnt[[i]], data = df_fin_15_rnt, family = gaussian)

  x <- t.test(df_fin_7_rnt[[i]], df_fin_17_rnt[[i]], var.equal = TRUE)
  p_val <- x[[3]]
  names(p_val) <- i
  t_test_p <- c(t_test_p, p_val)

}
t_test_p

#ANOVA of rank normalised values
aov_test_p <- c()
df_fin_rnt <- df_fin
for (i in mnames) {
  df_fin_rnt[[i]] <- as.numeric(df_fin_rnt[[i]])
  df_fin_rnt[[i]] <- rntransform(formula = df_fin_rnt[[i]], data = df_fin_rnt, family = gaussian)

  fom <- formula(paste("`", i, "`", " ~ ", "age", sep = ""))
  x <- aov(fom, data = df_fin_rnt)
  p_val <- summary(x)[[1]][[5]][1]
  names(p_val) <- i
  aov_test_p <- c(aov_test_p, p_val)
}

length(aov_test_p[aov_test_p < 0.05]) #217

####################################################
# Looking at the summary statistics of each metab ## - only using rank normalised data!
####################################################
#Function below takes the summary data from each metabolite and puts it into a new data frame
print_summ_stats <- function(df, columns) {
  summ_list <- list()
  for (i in columns) {
    met_name <- colnames(df[, i])
    summ <- summary(df[[i]])
    stats <- as.numeric(summ)
    names(stats) <- names(summ)
    summ_list[[met_name]] <- stats
  }
  df_2 <- do.call(rbind, summ_list)
  df_2 <- as.data.frame(df_2)
  df_2[["Metabolite"]] <- rownames(df_2)
  df_2
}

rnt = TRUE
rnt = FALSE

name <- ifelse(rnt, "_rnt", "")
if (rnt) {
  df <- df_fin_rnt
} else {
  df <- df_fin
}
#Create summary stats for the combined dataset and for data at each age
summ_stats_comb <- print_summ_stats(df, columns = 2:231)
summ_stat_list <- list()
for (i in unique(df[["age"]])) {
  df_temp <- subset(df, age == i)
  summ_stat_list[[i]] <- print_summ_stats(df_temp, columns = 2:231)
  colnames(summ_stat_list[[i]])[-7] <- paste(colnames(summ_stat_list[[i]])[-7], as.character(i), sep = "_")
}

#Make a dataset with summary stats at all ages in it
summ_stat_all <- left_join(summ_stat_list[[7]], summ_stat_list[[15]]) %>%
  left_join(., summ_stat_list[[17]]) %>%
  select(Metabolite, everything())

#Means only dataset
Means <- summ_stat_all %>%
  select(Metabolite, one_of(colnames(.)[grep("Mean", colnames(.))])) %>%
  mutate(rel_diff_7_17 = Mean_7 / Mean_17) %>%
  mutate(rel_diff_7_15 = Mean_7 / Mean_15)

#Medians only dataset
Medians <- summ_stat_all %>%
  select(Metabolite, one_of(colnames(.)[grep("Median", colnames(.))])) %>%
  mutate(rel_diff_7_17 = Median_7 / Median_17) %>%
  mutate(rel_diff_7_15 = Median_7 / Median_15)


#Print out the summ_stat_all table, put it into excel, highlight the crazy results manually
write.table(summ_stat_all, file = paste("summary_stats_at_all_ages", name, ".txt", sep = ""), quote = F, col.names = T, row.names = F, sep = "\t")


Means %>% 
  filter(Metabolite %in% nr_mnames) %>%
  filter(rel_diff_7_17 > 2 | rel_diff_7_17 < 0.5)

Means %>%
  filter(!(Metabolite %in% nr_mnames)) %>%
  filter(rel_diff_7_17 > 2 | rel_diff_7_17 < 0.5)

Means %>% 
  filter(Metabolite %in% nr_mnames) %>%
  filter(rel_diff_7_15 > 2 | rel_diff_7_15 < 0.5)

Means %>%
  filter(!(Metabolite %in% nr_mnames)) %>%
  filter(rel_diff_7_15 > 2 | rel_diff_7_15 < 0.5)



Medians %>%
  filter(Metabolite %in% nr_mnames) %>%
  filter(rel_diff_7_17 > 2 | rel_diff_7_17 < 0.5)

Medians %>%
  filter(!(Metabolite %in% nr_mnames)) %>%
  filter(rel_diff_7_17 > 2 | rel_diff_7_17 < 0.5)

Medians %>%
  filter(Metabolite %in% nr_mnames) %>%
  filter(rel_diff_7_15 > 2 | rel_diff_7_15 < 0.5)

Medians %>%
  filter(!(Metabolite %in% nr_mnames)) %>%
  filter(rel_diff_7_15 > 2 | rel_diff_7_15 < 0.5)



summ_stat_all %>%
  filter(Metabolite %in% c("l-VLDL-PL", "l-VLDL-L", "l-VLDL-PL_P"))

nr_means <- Means %>%
  filter(Metabolite %in% nr_mnames)

nr_medians <- Medians %>%
  filter(Metabolite %in% nr_mnames)

p <- ggplot(nr_means) +
  geom_point(aes(x = Mean_7, y = Mean_17)) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  theme(text = element_text(size = 20))
if (!rnt) {  
  p <- p + coord_cartesian(xlim = 0:2, ylim = 0:2)
}
pdf(paste("NR_Means_plot", name, ".pdf", sep = ""), width = 15, height = 10)
print(p)
dev.off()

p <- ggplot(nr_medians) +
  geom_point(aes(x = Median_7, y = Median_17)) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  theme(text = element_text(size = 20))

pdf(paste("NR_Medians_plot", name, ".pdf", sep = ""), width = 15, height = 10)
print(p)
dev.off()

ggplot(nr_medians) +
  geom_point(aes(x = Median_15, y = Median_17)) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  theme(text = element_text(size = 20))

ggplot(Means) +
  geom_point(aes(x = Mean_7, y = Mean_15)) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  coord_cartesian(xlim = 0:2, ylim = 0:2)

ggplot(Means) +
  geom_bar(aes(x = Metabolite, y = rel_diff_7_17), stat = "identity") +
  coord_cartesian(ylim = 0:3)



#######

summ_stat_list_2 <- list()
for (i in unique(df_fin[["age"]])) {
  df_temp <- subset(df_fin, age == i)
  summ_stat_list_2[[i]] <- print_summ_stats(df_temp, columns = 2:231)
  summ_stat_list_2[[i]][["age"]] <- i
}

summ_stat_all_2 <- do.call(rbind, summ_stat_list_2)
summ_stat_all_2

ggplot(summ_stat_all_2) +
  geom_point(aes(x = age, y = Mean)) +
  coord_cartesian(ylim = 0:2)


summ_stat_all <- do.call(cbind, summ_stat_list)
rownames(summ_stat_all) <- NULL
ylim1 <- boxplot.stats(summ_stat_all$Mean)$stats[c(1, 5)]


ggplot(summ_stat_all) +
  geom_boxplot(aes(x = as.factor(age), y = Mean)) +
  coord_cartesian(ylim = ylim1*1.05)

df_large <- df_f7
df_small <- test_df
col_ID <- 1

test_df2 <- sample_random_IDs(df_f7, test_df, 1, col_diff = "age")
length(unique(test_df2[[1]]))

table(test_df2[["age"]])


################################# 
# Checking main analysis ########
#################################

#rank normalise the metabolites

df_ages_list <- list()
for (i in unique(df_fin[["age"]])) {
  df_temp <- filter(df_fin, age == i)
  for (j in colnames(df_fin)[-c(1,232)]) {
    df_temp[[j]] <- as.numeric(df_temp[[j]])
    df_temp[[j]] <- rntransform(formula = df_temp[[j]], data = df_temp, family = gaussian)
  }
  df_ages_list[[i]] <- df_temp
}

df_fin <- do.call(rbind, df_ages_list)

mnames <- colnames(df_fin)[-c(1, 232)]

nr_mnames <- mnames[-grep("_P|/", mnames)]


#the sorted SNP file
datafile_SNPs <- read_dta("inputs/genotype_data.dta") #was New_Genotype_Sorted_dosage_genotype_data

#the SNPs and their weightings
SNP_Ws <- read_dta("inputs/SNP_weightings.dta") #was SNPs and log_odd_weightings_new_genotype

#########################################################
### 2. Produce SNP weightings  and sort out the SNPs ### 
#########################################################

##NB: Weightings are based on the log(OR) of each SNP with regards to CAD

#produce a new SNP file that can be altered and remove the IDs
#datafile_SNPs_W <- datafile_SNPs[-c(1,2)]
datafile_SNPs_W <- datafile_SNPs
datafile_SNPs_W <- as.data.frame(datafile_SNPs)
#str(datafile_SNPs_W)
#str(SNP_Ws)
rev_SNPs <- c("rs56062135", "rs8042271",  "rs180803", "rs11206510", "rs9970807",  "rs7528419",  "rs6689306",  "rs4593108",  "rs72689147", "rs6903956",  "rs17609940", "rs56336142", "rs12202017", "rs10953541", "rs11556924", "rs264",  "rs2954029",  "rs11191416", "rs2128739",  "rs964184", "rs3184504",  "rs12936587", "rs56289821")
#This loop changes the coding of the alleles so that the dosages all reflect the dosage of effect
#allele received by each individual then weights each SNP dosage
#To make this loop work need to change the tibbles into data frames
SNP_Ws <- as.data.frame(SNP_Ws)
for (i in colnames(datafile_SNPs_W)[-c(1, 2)]) {
  #Changes allele codings
  if (i %in% rev_SNPs) {
    datafile_SNPs_W[, i] <- 2 - datafile_SNPs_W[, i]
  }
  #Weights the SNPs
  weight <- subset(SNP_Ws, leadvariant == i)
  new_col_name <- paste(i, "_w", sep = "")
  datafile_SNPs_W[, new_col_name] <- datafile_SNPs_W[, i] * weight[, 2]
}
datafile_SNPs_W <- datafile_SNPs_W %>%
  mutate(u_ID = paste(cidB9999, qlet, sep = "_")) %>%
  dplyr::select(u_ID, everything())

####

d <- left_join(df_fin, datafile_SNPs_W) %>%
  dplyr::select(u_ID, cidB9999, qlet, age, everything())
#d2 <- d

#Drop the SNPs as decided before - including the indels - check this - will need to try and extract the indels!
bad_SNPs <- c("rs11830157_w", "rs12976411_w", "rs200200214_w", "rs199884676_w", "rs201810558_w", "rs201477372_w")

d <- dplyr::select(d, -one_of(bad_SNPs))

#d <- d[, -which(colnames(d) %in%  bad_SNPs)]

#extract SNP names - as well as lipid SNPs 
SNPs <-  colnames(d)[str_detect(colnames(d), "[\\d]_w$")]

#save the name of the SNP that has a log(OR) of 0
log_odd_SNP <- "rs6903956_w"

#lipid SNPs are those that fall within the gene regions of defined lipoprotein genes
lipid_SNPs <- c("rs11206510_w", "rs264_w", "rs1412444_w", "rs56289821_w", "rs4420638_w", "rs7528419_w", "rs964184_w")

#Produce a list of the SNPs that aren't directly related to lipoproteins 
non_lipid_SNPs <- SNPs[!(SNPs %in% lipid_SNPs)]

#order the data by unique identifiers
d <- arrange(d, cidB9999, qlet)

#########################################################
### 4. Generate the CAD and lipid risk scores #########
#########################################################
d[, "CAD_score"] <- rowSums(d[, SNPs])

d[, "lipid_score"] <- rowSums(d[, lipid_SNPs])

d[, "non_lipid_score"] <- rowSums(d[, non_lipid_SNPs])

##############################################
### 5. Do the linear regression analysis #####
##############################################

#source the linear regression function - mostly made by Qin Wang
source("R/Linear_regression_func.R")

#drop the SNP with a log odds ratio of 0
SNPs <- SNPs %>%
  .[!(. %in% log_odd_SNP)]

d <- select(d, -one_of(log_odd_SNP))

#Use lapply to do the linear regression analysis of each SNP on each metabolite - should produce a list of data frames

#Change d back into a data frame so it can be used with the linearRegress function
#indi_SNP_results <- list()
#indi_SNP_results_nr <- list()
CAD_score_lr_nr <- list()
lipid_score_lr_nr <- list()
non_lipid_score_lr_nr <- list()
for (i in unique(d[["age"]])) {
  d_temp <- filter(d, age == i)
  age <- as.character(i)
  #indi_SNP_results[[age]] <- lapply(SNPs, function(x) {linearRegress(x, mnames, d_temp)})
  #indi_SNP_results_nr[[age]] <- lapply(SNPs, function(x) {linearRegress(x, nr_mnames, d_temp)})
  CAD_score_lr_nr[[age]] <- linearRegress('CAD_score', nr_mnames, d_temp)
  lipid_score_lr_nr[[age]] <- linearRegress('lipid_score', nr_mnames, d_temp)
  non_lipid_score_lr_nr[[age]] <- linearRegress('non_lipid_score', nr_mnames, d_temp)
}

CAD_score_lr_nr_tot <- linearRegress('CAD_score', nr_mnames, d)
CAD_score_lr_nr_tot_age <- linearRegress('CAD_score', nr_mnames, d, covariate.name = "age")

############################################################################
#### Not sure what to do with this bit - look into it ####################
############################################################################

##For report
#x <- as.matrix(summary(CAD_score_lr_nr$Estimate))
#x <- cbind(x, summary(lipid_score_2_lr_nr$Estimate))
#x <- cbind(x, summary(non_lipid_score_lr_nr$Estimate))
#x <- cbind(x, summary(best_SNPs_score_lr_nr$Estimate))

#x <- t(x)

#rownames(x) <- c("CAD_score", "lipid_score", "non_lipid_score", "best_SNPs_score")

#x
#write.xlsx(x, "scores_vs_metab_coef_sum.xlsx")


#input the names of SNPs as the name of each data frame

SNP_names <- unlist(strsplit(SNPs, "_w"))

names(indi_SNP_results[[1]]) <- SNP_names

names(indi_SNP_results_nr[[1]]) <- SNP_names


############################################################################
###6. Extract the SNPs and metabolites with significant associations #######
############################################################################

indi_SNPs_p <- 0.05 / (length(SNPs) * length(mnames))
indi_SNPs_p_nr <- 0.05 / (length(SNPs) * length(nr_mnames))

#Add Bonferroni and FDR corrected p values to each of the results 
for(i in 1:length(indi_SNP_results)) {
  indi_SNP_results[[i]]$Bonferroni <- p.adjust(indi_SNP_results[[i]]$`Pr(>|t|)`, method = "bonferroni", n = (length(SNPs) * length(mnames)))
  indi_SNP_results[[i]]$FDR <- p.adjust(indi_SNP_results[[i]]$`Pr(>|t|)`, method = "fdr", n = (length(SNPs) * length(mnames)))
}

for(i in 1:length(indi_SNP_results_nr)) {
  indi_SNP_results_nr[[i]]$Bonferroni <- p.adjust(indi_SNP_results_nr[[i]]$`Pr(>|t|)`, method = "bonferroni", n = (length(SNPs) * length(nr_mnames)))
  indi_SNP_results_nr[[i]]$FDR <- p.adjust(indi_SNP_results_nr[[i]]$`Pr(>|t|)`, method = "fdr", n = (length(SNPs) * length(nr_mnames)))
}

p_score <- 0.05/length(mnames)
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

LP <- extract_sig_hits(indi_SNP_results)
LP_FDR <- extract_sig_hits(indi_SNP_results, type = "fdr")

LP_nr <- extract_sig_hits(indi_SNP_results_nr)
LP_nr_FDR <- extract_sig_hits(indi_SNP_results_nr, type = "fdr")

#############################################################################
names_LP <- names(LP)
names_LP_FDR <- names(LP_FDR)

names_LP_nr <- names(LP_nr)
names_LP_nr_FDR <- names(LP_nr_FDR)

#############################################################################
###
# JUST CAD SCORE
###
CAD_score_lr_nr_tot
CAD_score_lr_nr_tot_age
res_out <- cbind(CAD_score_lr_nr[["7"]], CAD_score_lr_nr[["15"]]$`Pr(>|t|)`, CAD_score_lr_nr[["17"]]$`Pr(>|t|)`, CAD_score_lr_nr_tot$`Pr(>|t|)`, CAD_score_lr_nr_tot_age$`Pr(>|t|)`)

  colnames(res_out)[8:11] <- c("P_val_15", "P_val_17", "P_val_tot", "P_val_cov_age")
  res_out <- res_out %>%
    mutate(exp_p = (1:nrow(.))/(nrow(.)+1)) %>%
    mutate(x = -log10(exp_p)) %>%
    mutate(y.7 = -log10(sort(`Pr(>|t|)`))) %>%
    mutate(y.15 = -log10(sort(P_val_15))) %>%
    mutate(y.17 = -log10(sort(P_val_17))) %>%
    mutate(y.tot = -log10(sort(P_val_tot))) %>%
    mutate(y.cov_age = -log10(sort(P_val_cov_age)))

  l <- list(estlambda(res_out$`Pr(>|t|)`), estlambda(res_out$P_val_15), estlambda(res_out$P_val_17), estlambda(res_out$P_val_tot), estlambda(res_out$P_val_cov_age))
  nom <- paste("CAD scores at ages 7, 15 and 17 vs. all metabolites", "\n", "number of tests = ", nrow(res_out), "\nAge 7: lambda = ", round(l[[1]]$estimate, 3), ", se = ", round(l[[1]]$se, 3), "\nAge 15: lambda = ", round(l[[2]]$estimate, 3), ", se = ", round(l[[2]]$se, 3), "\nAge 17: lambda = ", round(l[[3]]$estimate, 3), ", se = ", round(l[[3]]$se, 3), "\nAll_ages: lambda = ", round(l[[4]]$estimate, 3), ", se = ", round(l[[4]]$se, 3), "\nAge_covar: lambda = ", round(l[[5]]$estimate, 3), ", se = ", round(l[[5]]$se, 3), sep = "")
  p <- ggplot(res_out, aes(x = x, y = y.7)) +
    geom_point(aes(colour = "7")) +
    geom_point(aes(x = x, y = y.15, colour = "15")) +
    geom_point(aes(x = x, y = y.17, colour = "17")) +
    geom_point(aes(x = x, y = y.tot, colour = "all_ages")) +
    geom_point(aes(x = x, y = y.cov_age, colour = "cov_age")) +
    geom_abline(slope = 1, intercept = 0, colour = "red") +
    ggtitle(nom) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = expression(Expected ~ ~-log[10](P)), y = expression(Observed ~ ~-log[10](P)), colour = "Age")
pdf("outputs/CAD_GRS_age_comparison_3.pdf", width = 15, height = 10)
print(p)
dev.off()


########################################
### 9. Produce forest plots ############
########################################

source("R/forest_plot_function_Tom's.R")

################# No ratios metabolites forest plots ###########################################################

####COMPARING CAD_SCORES AND LIPID SCORES 
result <- list(CAD_score_lr_nr[["7"]], CAD_score_lr_nr[["15"]], CAD_score_lr_nr[["17"]])
names(result) <- c("7", "15", "17")
pdf("outputs/forests/CAD_score_age_comparison_nrmetabs_forest.pdf", width = 15, height = 10)
forestPlot(result, columns = 4)
dev.off()










