####################
# New GRS analysis #
####################

#Read in variants 
sig_SNPs <- read.delim("outputs/significant_SNPs.txt", stringsAsFactors = F, header = F)
#sig_SNPs <- read.delim("outputs/FDR_significant_SNPs.txt")
sig_SNPs <- paste(sig_SNPs[[1]], "_w", sep = "")
non_sig_SNPs <- SNPs[!(SNPs %in% sig_SNPs)]

#Produce a score based on SNPs that associated with 1 or more metabolite
d[["sig_score"]] <- rowSums(d[, sig_SNPs])
#Produce a score based on SNPs that didn't associate with any metabolites
d[["non_sig_score"]] <- rowSums(d[, non_sig_SNPs])

sig_score_lr_nr <- linearRegress('sig_score', nr_mnames, d)

non_sig_score_lr_nr <- linearRegress('non_sig_score', nr_mnames, d)

##########
# Plots ##
##########

##########
# Forest #
##########
#source("R/forest_plot_function_Tom's.R")
source("R/Forest_plot_functions.R")
non_sig_score_lr_nr <- arrange(non_sig_score_lr_nr, Metabolite)
sig_score_lr_nr <- arrange(sig_score_lr_nr, Metabolite)
CAD_score_lr_nr <- arrange(CAD_score_lr_nr, Metabolite)
stopifnot(nrow(CAD_score_lr_nr) == nrow(sig_score_lr_nr) && nrow(CAD_score_lr_nr) == nrow(sig_score_lr_nr))

facet_var <- facet_var_gen(CAD_score_lr_nr, col_num = 4, group_num = 3)
forest_dat <- rbind(CAD_score_lr_nr[, 1:7], sig_score_lr_nr, non_sig_score_lr_nr) %>%
  mutate(score = rep(c("CAD", "sig", "non_sig"), each = nrow(CAD_score_lr_nr))) %>%
  mutate(facet_var = facet_var)

forest_dat[["facet_var"]] <- as.factor(forest_dat[["facet_var"]])
pdf("outputs/forests/all_scores_nrmetabs_forest.pdf", width = 15, height = 10)
forest_plot(forest_dat, col_num = 4, group = "score", y_axis_var = "Metabolite", units = "Beta coefficient (95% CI)")
dev.off()

######
# QQ #
######

res_out <- cbind(CAD_score_lr_nr, sig_score_lr_nr$`Pr(>|t|)`, non_sig_score_lr_nr$`Pr(>|t|)`)
colnames(res_out)[10:11] <- c("sig_P", "non_sig_P")
res_out <- res_out %>%
  mutate(exp_p = (1:nrow(.))/(nrow(.)+1)) %>%
  mutate(x = -log10(exp_p)) %>%
  mutate(y.CAD = -log10(sort(`Pr(>|t|)`))) %>%
  mutate(y.sig = -log10(sort(sig_P))) %>%
  mutate(y.non_sig = -log10(sort(non_sig_P)))

l <- list(estlambda(res_out$`Pr(>|t|)`), estlambda(res_out$sig_P), estlambda(res_out$non_sig_P))
nom <- paste("CAD scores vs. all metabolites", "\n", "number of tests = ", nrow(res_out), "\nTotal score: lambda = ", round(l[[1]]$estimate, 3), ", se = ", round(l[[1]]$se, 3), "\nsig score: lambda = ", round(l[[2]]$estimate, 3), ", se = ", round(l[[2]]$se, 3), "\nNon-sig score: lambda = ", round(l[[3]]$estimate, 3), ", se = ", round(l[[3]]$se, 3), sep = "")
p <- ggplot(res_out, aes(x = x, y = y.CAD)) +
  geom_point(aes(colour = "Total GRS (56 variants)")) +
  geom_point(aes(x = x, y = y.sig, colour = "New GRS (7 variants)")) +
  geom_point(aes(x = x, y = y.non_sig, colour = "Residual GRS (49 variants)")) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  scale_colour_discrete(breaks = c("New GRS (7 variants)", "Total GRS (56 variants)", "Residual GRS (49 variants)")) +
  #ggtitle(nom) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 30), legend.text = element_text(size = 20)) +
  labs(x = expression(Expected ~ ~-log[10](P)), y = expression(Observed ~ ~-log[10](P)), colour = "Genetic risk score")
pdf("outputs/3_scores_vs_metabs_qq.pdf", width = 20, height = 10)
print(p)
dev.off()

##################
# For DOHaD talk #
##################
res_out <- cbind(CAD_score_lr_nr, sig_score_lr_nr$`Pr(>|t|)`, non_sig_score_lr_nr$`Pr(>|t|)`)
colnames(res_out)[10:11] <- c("sig_P", "non_sig_P")
res_out <- res_out %>%
  mutate(exp_p = (1:nrow(.))/(nrow(.)+1)) %>%
  mutate(x = -log10(exp_p)) %>%
  mutate(y.CAD = -log10(sort(`Pr(>|t|)`))) %>%
  mutate(y.sig = -log10(sort(sig_P))) %>%
  mutate(y.non_sig = -log10(sort(non_sig_P)))

l <- list(estlambda(res_out$`Pr(>|t|)`), estlambda(res_out$sig_P), estlambda(res_out$non_sig_P))
nom <- paste("CAD scores vs. all metabolites", "\n", "number of tests = ", nrow(res_out), "\nTotal score: lambda = ", round(l[[1]]$estimate, 3), ", se = ", round(l[[1]]$se, 3), "\nsig score: lambda = ", round(l[[2]]$estimate, 3), ", se = ", round(l[[2]]$se, 3), "\nNon-sig score: lambda = ", round(l[[3]]$estimate, 3), ", se = ", round(l[[3]]$se, 3), sep = "")
p <- ggplot(res_out, aes(x = x, y = y.CAD)) +
  geom_point(aes(colour = "Total GRS (56 variants)")) +
  #geom_point(aes(x = x, y = y.sig, colour = "New GRS (7 variants)")) +
  geom_point(aes(x = x, y = y.non_sig, colour = "Residual GRS (49 variants)")) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  scale_colour_discrete(breaks = c("New GRS (7 variants)", "Total GRS (56 variants)", "Residual GRS (49 variants)")) +
  #ggtitle(nom) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 30), legend.text = element_text(size = 20), axis.line = element_line(colour = "black")) +
  labs(x = expression(Expected ~ ~-log[10](P)), y = expression(Observed ~ ~-log[10](P)), colour = "Genetic risk score")
pdf("outputs/2_scores_vs_metabs_qq.pdf", width = 20, height = 10)
print(p)
dev.off()


print("Analysis complete")

#####################
# Association tests #
#####################
sig_score_assoc <- sig_score_lr_nr %>%
  mutate(Bonferroni = p.adjust(`Pr(>|t|)`, method = "bonferroni")) %>%
  mutate(FDR = p.adjust(`Pr(>|t|)`, method = "fdr")) %>%
  filter(FDR < 0.05)

summary(sig_score_lr_nr)
median(abs(sig_score_lr_nr$Estimate))

non_sig_score_assoc <- non_sig_score_lr_nr %>%
  mutate(Bonferroni = p.adjust(`Pr(>|t|)`, method = "bonferroni")) %>%
  mutate(FDR = p.adjust(`Pr(>|t|)`, method = "fdr")) %>%
  filter(FDR < 0.05)





