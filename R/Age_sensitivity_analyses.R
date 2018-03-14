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
forest_dat$age <- factor(forest_dat$age, levels = c(7, 15, 17))
pdf("outputs/forests/Age_diff_met_conc_no_particle_forest.pdf", width = 15, height = 10)
print(forest_plot(forest_dat, col_num = 4, group = "age", y_axis = "Metabolite", units = "Log10(median metabolite concentration)", null_at = 0))
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
forest_dat$age <- factor(forest_dat$age, levels = c(7, 15, 17))
pdf("outputs/forests/Age_diff_met_conc_particle_only_forest.pdf", width = 15, height = 10)
print(forest_plot(forest_dat, col_num = 2, group = "age", y_axis = "Metabolite", units = "Log10(median metabolite concentration)", null_at = 0))
dev.off()

# ----------------------------------------------------------------------------
# GRS-metabolite association 
# ----------------------------------------------------------------------------

full_dat <- do.call(rbind, age_diff)
full_dat$age <- c(rep(7, length(nr_mnames)), rep(15, length(nr_mnames)), rep(17, length(nr_mnames)))

kw_age_test <- kruskal.test(Estimate ~ age, data = full_dat) 
kw_age_test #p = 0.8227 therefore no evidence of difference between log10 median values of each metabolite at different age groups

colnames(met_dat)

if (!("CAD_score" %in% colnames(d))) {
	d[["CAD_score"]] <- rowSums(d[, SNPs])
}

#Comparison of metabolite-GRS associations
age_diff_reg <-  list()
for (i in unique(d[["age"]])) {
	temp_dat <- filter(d, age == i)
	age_diff_reg[[as.character(i)]] <- linearRegress('CAD_score', nr_mnames, temp_dat)
}
result <- age_diff_reg

stopifnot(result[[1]][["Metabolite"]] == result[[2]][["Metabolite"]] && result[[1]][["Metabolite"]] == result[[3]][["Metabolite"]])

facet_var <- facet_var_gen(result[[1]], col_num = 4, group_num = 3)

forest_dat <- do.call(rbind, result) %>%
	mutate(age = rep(c("17", "15", "7"), each = nrow(result[[1]]))) %>%
	mutate(facet_var = facet_var)

forest_dat[["facet_var"]] <- as.factor(forest_dat[["facet_var"]])
forest_dat$age <- factor(forest_dat$age, levels = c(7, 15, 17))
pdf("outputs/forests/Age_diff_GRS_met_forest.pdf", width = 15, height = 10)
print(forest_plot(forest_dat, col_num = 4, group = "age", y_axis = "Metabolite", units = "Beta coefficient (95% CI)", null_at = 0))
dev.off()

full_dat2 <- do.call(rbind, age_diff_reg)
full_dat2$age <- c(rep(17, 149), rep(15, 149), rep(7, 149))

kw_age_test2 <- kruskal.test(Estimate ~ age, data = full_dat2) 
kw_age_test2 #p = 0.007643 therefore evidence of difference between age groups 
dunn_age_test <- dunnTest(Estimate ~ as.factor(age), data = full_dat2, method = "bh")

# ----------------------------------------------------------------------------
# qq-plot for association between CAD-GRS and metabolites at each age
# ----------------------------------------------------------------------------
# full_dat2

# rownames(full_dat2) <- NULL
# # QQ 
# res_out <- select(full_dat2, age, `Pr(>|t|)`, Metabolite) %>%
# 	spread(key = age, value = `Pr(>|t|)`) %>%
# 	mutate(exp_p = (1:nrow(.))/(nrow(.)+1)) %>%
# 	mutate(x = -log10(exp_p)) %>%
# 	mutate(y.7 = -log10(sort(`7`))) %>%
# 	mutate(y.15 = -log10(sort(`15`))) %>%
# 	mutate(y.17 = -log10(sort(`17`)))

# #l <- list(estlambda(res_out$`Pr(>|t|)`), estlambda(res_out$sig_P), estlambda(res_out$non_sig_P))
# #nom <- paste("CAD scores vs. all metabolites", "\n", "number of tests = ", nrow(res_out), "\nTotal score: lambda = ", round(l[[1]]$estimate, 3), ", se = ", round(l[[1]]$se, 3), "\nsig score: lambda = ", round(l[[2]]$estimate, 3), ", se = ", round(l[[2]]$se, 3), "\nNon-sig score: lambda = ", round(l[[3]]$estimate, 3), ", se = ", round(l[[3]]$se, 3), sep = "")
# p <- ggplot(res_out, aes(x = x, y = y.7)) +
#   geom_point(aes(colour = "7")) +
#   geom_point(aes(x = x, y = y.15, colour = "15")) +
#   geom_point(aes(x = x, y = y.17, colour = "17")) +
#   geom_abline(slope = 1, intercept = 0, colour = "red") +
#   scale_colour_discrete(breaks = c("7", "15", "17")) +
#   #ggtitle(nom) +
#   theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 30), legend.text = element_text(size = 15)) +
#   labs(x = expression(Expected ~ ~-log[10](P)), y = expression(Observed ~ ~-log[10](P)), colour = "Age")
# pdf(paste0("outputs/other/3_ages_CAD-GRS_vs_metabs_qq.pdf"), width = 20, height = 10)
# print(p)
# dev.off()

# ----------------------------------------------------------------------------
# Age group differences using biological groupings
# ----------------------------------------------------------------------------

res_out <- left_join(full_dat2, subset_df) %>%
	filter(group != "NA")

# Would be nice if could change this so that the ages are in order...
p <- ggplot(arrange(res_out, age), aes(x = reorder(subset, abs(Estimate), FUN = median), y = abs(Estimate), fill = as.character(sort(age)))) +
  geom_boxplot() +
  labs(x = "Particle size class", y = "Relative effect estimate", legend = "Age") +
  theme(text = element_text(size = 25), axis.title.x = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_discrete(labels=c("Small_HDL" = "Small HDL", "V_Large_HDL" = "Very large HDL", "Large_VLDL" = "Large VLDL", "Large_HDL" = "Large HDL", "Remnant_particles" = "Remnant particles", "LDL" = "LDL")) +
  scale_fill_discrete(name = "Age")
pdf("outputs/other/age_strat_lipoprotein_subclass_effect_comparison.pdf",width = 15, height = 10)
print(p)
dev.off()

kw_est <- vector(mode = "list", length = length(unique(res_out$subset)))
names(kw_est) <- unique(res_out$subset)
dunn_est <- vector(mode = "list", length = length(unique(res_out$subset)))
names(dunn_est) <- unique(res_out$subset)
for (i in unique(res_out$subset)) {
	temp_dat <- filter(res_out, subset == i)
	kw_est[[i]] <- kruskal.test(Estimate ~ age, data = temp_dat) 
	dunn_est[[i]] <- dunnTest(Estimate ~ as.factor(age), data = temp_dat, method = "bh")
}

kw_age_test2 <- kruskal.test(Estimate ~ age, data = full_dat2) 
kw_age_test2 #p = 0.007643 therefore evidence of difference between age groups 
dunn_age_test <- dunnTest(Estimate ~ as.factor(age), data = full_dat2, method = "bh")



# ----------------------------------------------------------------------------
# Kettunen GWAS res
# ----------------------------------------------------------------------------
# library(TwoSampleMR)
# ao <- available_outcomes()
# ao[grep("Kettunen", ao$author), ]

# # Going to need to change the results in the main analysis so they represent an increase in metabolite conc for sd increase in score
# # Or an increase in metabolite SD for an increase in score - try both
# # How can this be done when we don't have individual level data?
# # - Can't produce a score from summary stats...
# # - Could look at just SNPs that are "significantly" associated with mets in our study
# # - Could potentially meta-analyse the results for these SNPs...
# # Need to remove normalisation of metabolites - check Kettunen methods - i.e. did he log metabs etc. 

# library(MRInstruments)
# data(metab_qtls)
# head(metab_qtls)

# # Making a new metabolite dataset with log-transformed metabs
# met_dat <- df_main_not_transformed %>%
# 	dplyr::select(u_ID, age, everything())

# for (i in mnames) {
#   df_main_not_transformed[[i]] <- as.numeric(df_main_not_transformed[[i]])
#   df_main_not_transformed[[i]] <- log(df_main_not_transformed[[i]])
# }

# log_d <- left_join(df_main_not_transformed, datafile_SNPs_W) %>%
#   dplyr::select(u_ID, cidB9999, qlet, age, everything())

# log_d <- arrange(log_d, cidB9999, qlet)

# log_d[["CAD_score"]] <- rowSums(d[, SNPs])

# # Regress the metabolites on the CAD genetic risk score with age as a covariate
# CAD_score_lr_nr2 <- linearRegress('CAD_score', nr_mnames, log_d, "age")


# met_dat[, -c(1,2)]
# colnames(met_dat)


# # Producing per sd increase values
# sd_score <- sd(d[["CAD_score"]], na.rm = TRUE)
# CAD_score_lr_nr[["Estimate_per_sd"]] <- CAD_score_lr_nr[["Estimate"]] * sd_score



# GWAS_met_labs <- unique(metab_qtls$phenotype)
# GWAS_met_labs %in% nr_mnames
# sort(nr_mnames)
# # Not in nr_mnames:
# mets_to_remove <- c("Bis.DB.ratio", "Bis.FA.ratio", "CH2.DB.ratio", "CH2.in.FA", "DB.in.FA", "FAw79S")
# GWAS_met_labs <- GWAS_met_labs[!(GWAS_met_labs %in% mets_to_remove)]
# GWAS_met_labs <- gsub("\\.", "-", GWAS_met_labs)
# GWAS_met_labs[!(tolower(GWAS_met_labs) %in% tolower(nr_mnames))] # Metabs in GWAS, but not in our data
# # TotPG and otPUFA needs to be renamed
# GWAS_met_labs[GWAS_met_labs %in% "otPUFA"] <- "PUFA"
# GWAS_met_labs[GWAS_met_labs %in% "TotPG"] <- "Tot-PG"

# mets_to_keep <- nr_mnames[tolower(nr_mnames) %in% tolower(GWAS_met_labs)]

# SNPs_to_test <- SNP_names[SNP_names %in% metab_qtls[["SNP"]] & SNP_names %in% names(sig_nr_FDR)]

# SNP_sds <- vector(mode = "numeric", length = length(SNPs_to_test))
# names(SNP_sds) <- SNPs_to_test
# for (i in SNPs_to_test) {
# 	SNP_sds[i] <- sd(d[[i]], na.rm = TRUE)
# }

# # Test
# child_res <- indi_SNP_results_nr[[SNPs_to_test]] %>%
# 	filter(Metabolite %in% mets_to_keep) %>%
# 	mutate(Estimate_per_sd = Estimate * SNP_sds)






