# ----------------------------------------------------------------------------
# Age sensitivity analyses
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Setup of data
# ----------------------------------------------------------------------------

# Comparison of metabolite concentrations 
met_dat <- df_main_not_transformed %>%
	dplyr::select(u_ID, age, everything())
temp <- summary(met_dat[, 3])

data <- met_dat[, 3]

# Split the metabolite data into different age groups
age_diff <- list()
age_diff_not_log <- list()
for (i in unique(met_dat[["age"]])) {
	temp_dat <- filter(met_dat, age == i)
	medians <- sapply(dplyr::select(temp_dat, one_of(nr_mnames)), function(x) {median(x)})
	log_med <- log10(medians)
	low_quart <- sapply(dplyr::select(temp_dat, one_of(nr_mnames)), function(x) {as.numeric(quantile(x, na.rm = TRUE)[2])})
	log_low_quart <- log10(low_quart)
	up_quart <- sapply(dplyr::select(temp_dat, one_of(nr_mnames)), function(x) {as.numeric(quantile(x, na.rm = TRUE)[4])})
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


full_dat <- do.call(rbind, age_diff)
full_dat$age <- c(rep(7, length(nr_mnames)), rep(15, length(nr_mnames)), rep(17, length(nr_mnames)))

kw_age_test <- kruskal.test(Estimate ~ age, data = full_dat) 
kw_age_test #p = 0.8232 therefore no evidence of difference between log10 median values of each metabolite at different age groups

colnames(met_dat)

source("R/Forest_plot_functions.R")

# ----------------------------------------------------------------------------
# GRS-metabolite association 
# ----------------------------------------------------------------------------

if (!("CAD_score" %in% colnames(d))) {
	d[["CAD_score"]] <- rowSums(d[, SNPs])
}

#Comparison of metabolite-GRS associations
age_diff_reg <-  list()
for (i in unique(d[["tp"]])) {
	temp_dat <- filter(d, tp == i)
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
p <- forest_plot(forest_dat, col_num = 4, group = "age", y_axis = "Metabolite", units = "Beta coefficient (95% CI)", null_at = 0)
ggsave("outputs/forests/Age_diff_GRS_met_forest.pdf", plot = p, width = 15, height = 10, units = "in")

full_dat2 <- do.call(rbind, age_diff_reg)
full_dat2$age <- c(rep(17, length(nr_mnames)), rep(15, length(nr_mnames)), rep(7, length(nr_mnames)))

kw_age_test2 <- kruskal.test(Estimate ~ age, data = full_dat2) 
print(kw_age_test2) #p = 0.02417 therefore evidence of difference between age groups 
dunn_age_test <- dunnTest(Estimate ~ as.factor(age), data = full_dat2, method = "bh")

# ----------------------------------------------------------------------------
# Age group differences using biological groupings
# ----------------------------------------------------------------------------

res_out <- left_join(full_dat2, subset_df) %>%
	filter(group != "NA")
res_out$age <- factor(res_out$age, levels = c(7,15,17))

p <- ggplot(res_out, aes(x = reorder(subset, abs(Estimate), FUN = median), y = abs(Estimate), fill = age)) +
  geom_boxplot() +
  labs(x = "Particle size class", y = "Effect estimate", legend = "Age") +
  theme(text = element_text(size = 25), axis.title.x = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size = 20)) +
  scale_x_discrete(labels=c("Small_HDL" = "Small HDL", "Very_large_HDL" = "Very large HDL", "Large_VLDL" = "Large VLDL", "Large_HDL" = "Large HDL", "Atherogenic_non_LDL" = "Atherogenic non-LDL", "LDL" = "LDL")) +
  scale_fill_discrete(name = "Age", breaks = levels(res_out$age))

ggsave("outputs/other/age_strat_lipoprotein_subclass_effect_comparison.pdf", plot = p, width = 15, height = 10, units = "in")

# Kruskal-wallis - differences at each age
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


