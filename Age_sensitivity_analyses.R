#############################
# Age sensitivity analyses ##
#############################

#Make sure to complete all the parts in the masterfile before starting 

#So far it seems as though there is potentially a difference between those at aged 7, 15 and 17
# -little difference between 15 and 17
# -little difference between the medians of log10 transformed metabolite concs 
# -but difference between the associations with GRS 
# -tried to find out which metabolite-GRS associations differ by age, but can't seem to do a loop to find out...
# -graphs also seem to suggest subtle differences
# -SOLUTION: there are many people aged 7, could just have the main analysis done at that age and then subsequent analyses done at age 15 and 17




#Comparison of metabolite concentrations 
met_dat <- df_main_not_transformed %>%
	dplyr::select(u_ID, age, everything())
temp <- summary(met_dat[,3])
summary(lm(met_dat[,3] ~ met_dat[,4]))
str_split(temp[1, ], ":")[[1]][2]

data <- met_dat[, 3]


age_diff <- list()
for (i in unique(met_dat[["age"]])) {
	temp_dat <- filter(met_dat, age == i)
	medians <- sapply(select(temp_dat, one_of(nr_mnames)), function(x) {median(x)})
	log_med <- log10(medians)
	low_quart <- sapply(select(temp_dat, one_of(nr_mnames)), function(x) {as.numeric(quantile(x, na.rm = TRUE)[2])})
	log_low_quart <- log10(low_quart)
	up_quart <- sapply(select(temp_dat, one_of(nr_mnames)), function(x) {as.numeric(quantile(x, na.rm = TRUE)[4])})
	log_up_quart <- log10(up_quart)

	temp_out <- data.frame(log_med, log_low_quart, log_up_quart) %>%
		rownames_to_column()
	colnames(temp_out) <- c("Metabolite", "Estimate", "2.5 %", "97.5 %")
	temp_out <- filter(temp_out, !(Metabolite %in% c("glucose", "insulin"))) %>%
		mutate(`Pr(>|t|)` = NA)
	age_diff[[as.character(i)]] <- temp_out
}

particle_mets <- nr_mnames[grep("-P$", nr_mnames)]
result <- age_diff
for (i in names(result)) {
	result[[i]] <- result[[i]] %>%
		filter(!(Metabolite %in% particle_mets))
}
###PRODUCES A FOREST THAT ALLOWS VISUALISATION OF DATA - the metabolite conc differences are too large
# Still not completely right
pdf("outputs/Age_diff_met_conc_no_particle_forest.pdf", width = 15, height = 10)
forestPlot(result, columns = 4)
dev.off()

#--------- Particle sizes only

result <- age_diff
for (i in names(result)) {
	result[[i]] <- result[[i]] %>%
		filter(Metabolite %in% particle_mets)
}
###PRODUCES A FOREST THAT ALLOWS VISUALISATION OF DATA - the metabolite conc differences are too large
# Still not completely right
pdf("outputs/Age_diff_met_conc_particle_only_forest.pdf", width = 15, height = 10)
forestPlot(result, columns = 1)
dev.off()

#--------------------------------

full_dat <- do.call(rbind, age_diff)
full_dat$age <- c(rep(7, 149), rep(15, 149), rep(17, 149))

kw_age_test <- kruskal.test(Estimate ~ age, data = full_dat) 
kw_age_test #p = 0.8227 therefore no evidence of difference between age groups 

#Comparison of metabolite-GRS associations
age_diff_reg <-  list()
for (i in unique(d[["age"]])) {
	temp_dat <- filter(d, age == i) %>%
		dplyr::select(-insulin, -glucose)
	age_diff_reg[[as.character(i)]] <- linearRegress('CAD_score', nr_mnames[-c(150, 151)], temp_dat)
}
result <- age_diff_reg
pdf("outputs/Age_diff_GRS_met_forest.pdf", width = 15, height = 10)
forestPlot(result, columns = 4)
dev.off()

full_dat2 <- do.call(rbind, age_diff_reg)
full_dat2$age <- c(rep(17, 149), rep(15, 149), rep(7, 149))

kw_age_test2 <- kruskal.test(Estimate ~ age, data = full_dat2) 
aov_test <- aov(Estimate ~ as.factor(age), data = full_dat2)
summary(aov_test)
kw_age_test2 #p = 0.007643 therefore evidence of difference between age groups 
dunn_age_test <- dunnTest(Estimate ~ as.factor(age), data = full_dat2, method = "bh")

#For some reason the p vals are all the same - no residuals? 
met_dunn_tests <- list()
for (i in nr_mnames[!nr_mnames %in% c("insulin", "glucose")]) {
	print(i)
	temp_dat <- subset(full_dat2, Metabolite == i)
	temp_kw <- kruskal.test(Estimate ~ age, data = temp_dat)
	if (temp_kw[[3]] < 0.05) {
		print(paste("There is evidence to suggest the association between the GRS and ", i, " differs depending on age", sep = ""))
		temp_dunn <- dunnTest(Estimate ~ as.factor(age), data = temp_dat, method = "bh")
		met_dunn_tests[[i]] <- temp_dunn
	}
}
hist(age_diff_reg[[1]][["Estimate"]])
hist(age_diff_reg[[2]][["Estimate"]])
hist(age_diff_reg[[3]][["Estimate"]])










