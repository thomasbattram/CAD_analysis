# ----------------------------------------------------------------------------
# Age sensitivity analyses
# ----------------------------------------------------------------------------

# Make sure to complete all the parts in the masterfile before starting 

#So far it seems as though there is potentially a difference between those at aged 7, 15 and 17
# -little difference between 15 and 17
# -little difference between the medians of log10 transformed metabolite concs 
# -but difference between the associations with GRS 
# -tried to find out which metabolite-GRS associations differ by age, but can't seem to do a loop to find out...
# -graphs also seem to suggest subtle differences
# -SOLUTION: there are many people aged 7, could just have the main analysis done at that age and then subsequent analyses done at age 15 and 17

# ----------------------------------------------------------------------------
# Setup of data
# ----------------------------------------------------------------------------

# Comparison of metabolite concentrations 
met_dat <- df_main_not_transformed %>%
	dplyr::select(u_ID, age, everything())
temp <- summary(met_dat[, 3])
str_split(temp[1, ], ":")[[1]][2]

data <- met_dat[, 3]

# Split the metabolite data into different age groups
age_diff <- list()
age_diff_not_log <- list()
for (i in unique(met_dat[["age"]])) {
	temp_dat <- filter(met_dat, age == i)
	medians <- sapply(select(temp_dat, one_of(nr_mnames)), function(x) {median(x)})
	log_med <- log10(medians)
	low_quart <- sapply(select(temp_dat, one_of(nr_mnames)), function(x) {as.numeric(quantile(x, na.rm = TRUE)[2])})
	log_low_quart <- log10(low_quart)
	up_quart <- sapply(select(temp_dat, one_of(nr_mnames)), function(x) {as.numeric(quantile(x, na.rm = TRUE)[4])})
	log_up_quart <- log10(up_quart)

	temp_out_log <- data.frame(log_med, log_low_quart, log_up_quart) %>%
		rownames_to_column()
	colnames(temp_out_log) <- c("Metabolite", "Estimate", "2.5 %", "97.5 %")
	temp_out_log <- filter(temp_out_log, !(Metabolite %in% c("glucose", "insulin"))) %>%
		mutate(`Pr(>|t|)` = NA)
	age_diff[[as.character(i)]] <- temp_out_log
	
	temp_out <- data.frame(medians, low_quart, up_quart) %>%
		rownames_to_column()
	colnames(temp_out) <- c("Metabolite", "Estimate", "2.5 %", "97.5 %")
	temp_out <- filter(temp_out, !(Metabolite %in% c("glucose", "insulin"))) %>%
		mutate(`Pr(>|t|)` = NA)

	age_diff_not_log[[as.character(i)]] <- temp_out
}

# ----------------------------------------------------------------------------
# Metabolite concentration 
# ----------------------------------------------------------------------------

# Remove metabolites that are particle sizes - values too different to put on same graph
particle_mets <- nr_mnames[grep("-P$", nr_mnames)]
result <- age_diff
for (i in names(result)) {
	result[[i]] <- result[[i]] %>%
		filter(!(Metabolite %in% particle_mets))
}

stopifnot(result[[1]][["Metabolite"]] == result[[2]][["Metabolite"]] && result[[1]][["Metabolite"]] == result[[3]][["Metabolite"]])

source("R/Forest_plot_functions.R")

facet_var <- facet_var_gen(result[[1]], col_num = 4, group_num = 3)

forest_dat <- do.call(rbind, result) %>%
	mutate(age = rep(c("7", "15", "17"), each = nrow(result[[1]]))) %>%
	mutate(facet_var = facet_var)

# Could be changed still? 
forest_dat[["facet_var"]] <- as.factor(forest_dat[["facet_var"]])
pdf("outputs/forests/Age_diff_met_conc_no_particle_forest.pdf", width = 15, height = 10)
forest_plot(forest_dat, col_num = 4, group = "age", y_axis_var = "Metabolite", units = "Log10(median metabolite concentration)")
dev.off()

# Particle sizes only
result <- age_diff
for (i in names(result)) {
	result[[i]] <- result[[i]] %>%
		filter(Metabolite %in% particle_mets)
}
stopifnot(result[[1]][["Metabolite"]] == result[[2]][["Metabolite"]] && result[[1]][["Metabolite"]] == result[[3]][["Metabolite"]])

#Need to change facet_var_gen()
facet_var <- facet_var_gen(result[[1]], col_num = 2, group_num = 3)

forest_dat <- do.call(rbind, result) %>%
	mutate(age = rep(c("7", "15", "17"), each = nrow(result[[1]]))) %>%
	mutate(facet_var = facet_var)

# Could be changed still? 
forest_dat[["facet_var"]] <- as.factor(forest_dat[["facet_var"]])
pdf("outputs/forests/Age_diff_met_conc_particle_only_forest.pdf", width = 15, height = 10)
forest_plot(forest_dat, col_num = 2, group = "age", y_axis_var = "Metabolite", units = "Log10(median metabolite concentration)")
dev.off()

# ----------------------------------------------------------------------------
# GRS-metabolite association 
# ----------------------------------------------------------------------------

full_dat <- do.call(rbind, age_diff)
full_dat$age <- c(rep(7, 149), rep(15, 149), rep(17, 149))

kw_age_test <- kruskal.test(Estimate ~ age, data = full_dat) 
kw_age_test #p = 0.8227 therefore no evidence of difference between log10 median values of each metabolite at different age groups

colnames(met_dat)

if (!("CAD_score" %in% colnames(d))) {
	d[["CAD_score"]] <- rowSums(d[, SNPs])
}

#Comparison of metabolite-GRS associations
age_diff_reg <-  list()
for (i in unique(d[["age"]])) {
	temp_dat <- filter(d, age == i) %>%
		dplyr::select(-insulin, -glucose)
	age_diff_reg[[as.character(i)]] <- linearRegress('CAD_score', nr_mnames[-c(150, 151)], temp_dat)
}
result <- age_diff_reg

stopifnot(result[[1]][["Metabolite"]] == result[[2]][["Metabolite"]] && result[[1]][["Metabolite"]] == result[[3]][["Metabolite"]])

facet_var <- facet_var_gen(result[[1]], col_num = 4, group_num = 3)

forest_dat <- do.call(rbind, result) %>%
	mutate(age = rep(c("17", "15", "7"), each = nrow(result[[1]]))) %>%
	mutate(facet_var = facet_var)

forest_dat[["facet_var"]] <- as.factor(forest_dat[["facet_var"]])

pdf("outputs/forests/Age_diff_GRS_met_forest.pdf", width = 15, height = 10)
forest_plot(forest_dat, col_num = 4, group = "age", y_axis_var = "Metabolite", units = "Beta coefficient (95% CI)")
dev.off()

full_dat2 <- do.call(rbind, age_diff_reg)
full_dat2$age <- c(rep(17, 149), rep(15, 149), rep(7, 149))

kw_age_test2 <- kruskal.test(Estimate ~ age, data = full_dat2) 
kw_age_test2 #p = 0.007643 therefore evidence of difference between age groups 
dunn_age_test <- dunnTest(Estimate ~ as.factor(age), data = full_dat2, method = "bh")

# Kettunen GWAS res
library(TwoSampleMR)
ao <- available_outcomes()
ao[grep("Kettunen", ao$author), ]

# Going to need to change the results in the main analysis so they represent an increase in metabolite conc for sd increase in score
# Or an increase in metabolite SD for an increase in score - try both
# How can this be done when we don't have individual level data?
# - Can't produce a score from summary stats...
# - Could look at just SNPs that are "significantly" associated with mets in our study
# - Could potentially meta-analyse the results for these SNPs...
# Need to remove normalisation of metabolites - check Kettunen methods - i.e. did he log metabs etc. 

library(MRInstruments)
data(metab_qtls)
head(metab_qtls)

# Producing per sd increase values
sd_score <- sd(d[["CAD_score"]], na.rm = TRUE)
CAD_score_lr_nr[["Estimate_per_sd"]] <- CAD_score_lr_nr[["Estimate"]] * sd_score



GWAS_met_labs <- unique(metab_qtls$phenotype)
GWAS_met_labs %in% nr_mnames
sort(nr_mnames)
# Not in nr_mnames:
mets_to_remove <- c("Bis.DB.ratio", "Bis.FA.ratio", "CH2.DB.ratio", "CH2.in.FA", "DB.in.FA", "FAw79S")
GWAS_met_labs <- GWAS_met_labs[!(GWAS_met_labs %in% mets_to_remove)]
GWAS_met_labs <- gsub("\\.", "-", GWAS_met_labs)
GWAS_met_labs[!(tolower(GWAS_met_labs) %in% tolower(nr_mnames))] # Metabs in GWAS, but not in our data
# TotPG and otPUFA needs to be renamed
GWAS_met_labs[GWAS_met_labs %in% "otPUFA"] <- "PUFA"
GWAS_met_labs[GWAS_met_labs %in% "TotPG"] <- "Tot-PG"

mets_to_keep <- nr_mnames[tolower(nr_mnames) %in% tolower(GWAS_met_labs)]

SNPs_to_test <- SNP_names[SNP_names %in% metab_qtls[["SNP"]] & SNP_names %in% names(sig_nr_FDR)]

SNP_sds <- vector(mode = "numeric", length = length(SNPs_to_test))
names(SNP_sds) <- SNPs_to_test
for (i in SNPs_to_test) {
	SNP_sds[i] <- sd(d[[i]], na.rm = TRUE)
}

# Test
child_res <- indi_SNP_results_nr[[SNPs_to_test]] %>%
	filter(Metabolite %in% mets_to_keep) %>%
	mutate(Estimate_per_sd = Estimate * SNP_sds)






