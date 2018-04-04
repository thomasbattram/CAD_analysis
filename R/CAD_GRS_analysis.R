# ------------------------------------------------------------------
# CAD-GRS metabolite analysis
# ------------------------------------------------------------------

d[["CAD_score"]] <- rowSums(d[, SNPs])

# ------------------------------------------------------------------
# Do the linear regression analysis 
# ------------------------------------------------------------------

# Regress the metabolites on the CAD genetic risk score with age as a covariate
CAD_score_lr_nr <- linearRegress('CAD_score', nr_mnames, d, "age")

# ------------------------------------------------------------------
# CAD-GRS vs. metabolites qq 
# ------------------------------------------------------------------
# Create expected and observed p value columns (x and y)
res_out <- CAD_score_lr_nr %>%
  mutate(exp_p = (1:nrow(.))/(nrow(.)+1)) %>%
  mutate(x = -log10(exp_p)) %>%
  mutate(y = -log10(sort(`Pr(>|t|)`)))

l <- list(estlambda(res_out$`Pr(>|t|)`))
nom <- paste("CAD scores vs. all metabolites", "\n", "number of tests = ", nrow(res_out), "\nlambda = ", round(l[[1]]$estimate, 3), ", se = ", round(l[[1]]$se, 3), sep = "") 
p <- ggplot(res_out, aes(x = x, y = y)) +
  geom_point(colour = "blue") +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  #ggtitle(nom) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 30)) +
  labs(x = expression(Expected ~ ~-log[10](P)), y = expression(Observed ~ ~-log[10](P)))
pdf(paste0("outputs/other/", as.character(age), "/CAD_score_qq.pdf"), width = 15, height = 10)
print(p)
dev.off()

# ------------------------------------------------------------------
# CAD-GRS vs. grouped lipoproteins qq 
# ------------------------------------------------------------------

CAD_score_lr_nr <- left_join(CAD_score_lr_nr, subset_df)

CAD_score_lr_nr$group <- as.factor(CAD_score_lr_nr$group)
CAD_score_lr_nr$subset <- as.factor(CAD_score_lr_nr$subset)

# Remove metabolites not grouped (non-lipoproteins)
dat <- CAD_score_lr_nr[-which(is.na(CAD_score_lr_nr$group)), ]
str(dat)
dat$subset <- factor(dat$subset, levels = c("LDL", "Remnant_particles", "Large_VLDL", "Small_HDL", "Large_HDL", "V_Large_HDL"))

# Create expected and observed p value columns (x and y)
res_out <- dat %>%
  arrange(`Pr(>|t|)`) %>%
  mutate(exp_p = (1:nrow(.))/(nrow(.)+1)) %>%
  mutate(x = -log10(exp_p)) %>%
  mutate(y = -log10(`Pr(>|t|)`))

# Create the plot
l <- estlambda(res_out$`Pr(>|t|)`)
nom <- paste("CAD score vs. lipoprotein metabolites", "\n", "number of tests = ", nrow(res_out), "\nlambda = ", round(l$estimate, 3), ", se = ", round(l$se, 3), sep = "")
  p <- ggplot(res_out, aes(x = x, y = y, group = group, colour = subset)) +
    geom_point(aes(shape = factor(group)), size = 3) +
    geom_abline(slope = 1, intercept = 0, colour = "red") +
    #ggtitle(nom) +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 30), legend.text = element_text(size = 20), legend.title = element_blank()) +
    labs(x = expression(Expected ~ ~-log[10](P)), y = expression(Observed ~ ~-log[10](P))) +
    scale_colour_hue(breaks = levels(res_out[["subset"]])) #+
    #theme(legend.position = "right") +
    #scale_colour_manual(name = element_blank(),
    ###labels need to be in alphabetical order!!
    #labels = c("l-HDL", "l-VLDL", "LDL", "Rem P", "s-HDL", "vl-HDL"),
    #values = colours))

pdf(paste0("outputs/other/", as.character(age), "/CAD_score_vs_lipoprotein_subclasses_QQ.pdf"), width = 15, height = 10)
print(p)
dev.off()

# ------------------------------------------------------------------
# Differences between subsets 
# ------------------------------------------------------------------
# Boxplot looking at effect size differences between groups
pdf(paste0("outputs/other/", as.character(age), "/lipoprotein_subclass_effect_comparison.pdf"), width = 15, height = 10)
p <- ggplot(res_out) +
  geom_boxplot(aes(x = reorder(subset, abs(Estimate), FUN = median), y = abs(Estimate), fill = subset)) +
  labs(x = "Particle size class", y = "Effect estimate") +
  theme(text = element_text(size = 30), legend.position = "none", axis.title.x = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_discrete(labels=c("Small_HDL" = "Small HDL", "V_Large_HDL" = "Very large HDL", "Large_VLDL" = "Large VLDL", "Large_HDL" = "Large HDL", "Remnant_particles" = "Remnant particles", "LDL" = "LDL"))
print(p)
dev.off()

# Testing for normality within group estimates
shap_test <- vector(mode = "list", length = length(unique(res_out$subset)))
names(shap_test) <- unique(res_out$subset)
for (i in res_out$subset) {
  temp_dat <- filter(res_out, subset == i)
  shap_test[[i]] <- shapiro.test(temp_dat$Estimate)
}
shap_test
# The estimates within many subsets are not normally distributed - use kruskal-wallis test

kw_test <- kruskal.test(Estimate ~ subset, data = res_out) 
dunn_test <- dunnTest(Estimate ~ subset, data = res_out, method = "bh")
dunn_res <- dunn_test[[2]]
write.table(dunn_res, paste0("outputs/tables/", as.character(age), "/lipoprotein_subclass_kw_test.txt"), quote = F, col.names = T, row.names = F, sep = "\t")

# ------------------------------------------------------------------
# Association tables 
# ------------------------------------------------------------------
CAD_assoc <- CAD_score_lr_nr %>%
  mutate(Bonferroni = p.adjust(`Pr(>|t|)`, method = "bonferroni")) %>%
  mutate(FDR = p.adjust(`Pr(>|t|)`, method = "fdr")) %>%
  filter(FDR < 0.05)

length(lipoproteins)  
sum(lipoproteins %in% CAD_assoc$Metabolite)


CAD_tab <- arrange(CAD_score_lr_nr, `Pr(>|t|)`) %>%
  mutate(P = make_pretty(`Pr(>|t|)`, 3)) %>%
  mutate(FDR = make_pretty(p.adjust(`Pr(>|t|)`, method = "fdr"), 3)) %>%
  mutate(low_CI = make_pretty(`2.5 %`, 3)) %>%
  mutate(up_CI = make_pretty(`97.5 %`, 3)) %>%
  mutate(`Estimate (95% CI)` = paste0(make_pretty(Estimate, 3), " (", low_CI, ", ", up_CI, ")")) %>%
  dplyr::select(Metabolite, `Estimate (95% CI)`, P, FDR)

write.table(CAD_tab, file = paste0("outputs/tables/", as.character(age), "/CAD_GRS_metabolite_full_assoc.txt"), quote = F, col.names = T, row.names = F, sep = "\t")

median(abs(CAD_score_lr_nr$Estimate))

# ------------------------------------------------------------------
# Independent features
# ------------------------------------------------------------------

assoc_met <- filter(CAD_tab, FDR < 0.05) %>%
  dplyr::select(Metabolite)
assoc_met <- assoc_met[, 1]

assoc_tab <- dplyr::select(d, one_of(assoc_met))

D <-  as.dist(1 - abs(cor(assoc_tab)))
PRhoTree <- hclust(D, method = "complete")
k <- cutree(PRhoTree, h = 0.2) ## distances are in 1 - Pearson Rho distances
table(k) ## 20 Clusters

print("CAD-GRS analysis complete - please proceed to the Individual variant analysis")




