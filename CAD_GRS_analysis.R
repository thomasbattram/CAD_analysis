###############################
# CAD-GRS metabolite analysis #
###############################
#install.packages("FSA")
library(FSA)


d[["CAD_score"]] <- rowSums(d[, SNPs])

##############################################
### Do the linear regression analysis #####
##############################################

#source the linear regression function - mostly made by Qin Wang
source("R/Linear_regression_func.R")

#drop the SNP with a log odds ratio of 0
SNPs <- SNPs %>%
  .[!(. %in% log_odd_SNP)]

d <- select(d, -one_of(log_odd_SNP))

CAD_score_lr_nr <- linearRegress('CAD_score', nr_mnames, d, "age")

##############################
# CAD-GRS vs. metabolites qq #
##############################
res_out <- CAD_score_lr_nr
res_out <- res_out %>%
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
pdf("outputs/CAD_score_qq.pdf", width = 15, height = 10)
print(p)
dev.off()

#######################################
# CAD-GRS vs. grouped lipoproteins qq #
#######################################
#Generate the lipoprotein groups
lipoproteins <- c(nr_mnames[grep("-V?[HIL]DL", nr_mnames)], nr_mnames[grep("IDL", nr_mnames)])
est_tot <- lipoproteins[grep(pattern = "-C", lipoproteins)]
free <- lipoproteins[grep(pattern = "-FC", lipoproteins)]
cholesterol <- lipoproteins[lipoproteins %in% c(free, est_tot)]
phospholipid <- lipoproteins[grep(pattern = "-PL", lipoproteins)]
triglyceride <- lipoproteins[grep(pattern = "-TG", lipoproteins)]
tot_lipid <- lipoproteins[grep(pattern = "L-L", lipoproteins)]
particle <- lipoproteins[grep(pattern = "-P", lipoproteins)]
particle <- particle[-which(particle %in% phospholipid)]

groups <- list(cholesterol, phospholipid, triglyceride, tot_lipid, particle)
names(groups) <- c("cholesterol", "phospholipid", "triglyceride", "tot_lipid", "particle")

#Generate the subsets 
Large_VLDL <- c("xxl-VLDL", "xl-VLDL", "l-VLDL")
Remnant_particles <- c("m-VLDL", "s-VLDL", "xs-VLDL", "IDL")
LDL <- c("l-LDL", "m-LDL", "s-LDL")
V_Large_HDL <- c("xl-HDL") 
Large_HDL <- c("l-HDL", "m-HDL")
Small_HDL <- "s-HDL"

subsets <- list(Large_VLDL, Remnant_particles, LDL, V_Large_HDL, Large_HDL, Small_HDL)
names(subsets) <- c("Large_VLDL", "Remnant_particles", "LDL", "V_Large_HDL", "Large_HDL", "Small_HDL")

#Add groups and subsets to the data
for (i in names(groups)) {
  CAD_score_lr_nr[CAD_score_lr_nr[["Metabolite"]] %in% groups[[i]], "group"] <- i
}
#CAD_score_lr_nr$subset <- NA
groups_2 <- list()
for (j in names(subsets)) {
    n <- subsets[[j]]
    o <- list()
    for (m in n) {
      o[m] <- list(CAD_score_lr_nr[["Metabolite"]][grep(paste("\\b",m,sep=""), CAD_score_lr_nr[["Metabolite"]])])
    }
    names(o) <- NULL
    groups_2[j] <- list(unlist(o))
    for (k in names(groups_2)) {
      CAD_score_lr_nr[CAD_score_lr_nr[["Metabolite"]] %in% groups_2[[k]], "subset"] <- k
    }
  }
CAD_score_lr_nr$group <- as.factor(CAD_score_lr_nr$group)
CAD_score_lr_nr$subset <- as.factor(CAD_score_lr_nr$subset)

#Remove metabolites not grouped (non-lipoproteins)
dat <- CAD_score_lr_nr[-which(is.na(CAD_score_lr_nr$group)), ]

#Create expected and observed p value columns (x and y)
res_out <- dat %>%
  arrange(`Pr(>|t|)`) %>%
  mutate(exp_p = (1:nrow(.))/(nrow(.)+1)) %>%
  mutate(x = -log10(exp_p)) %>%
  mutate(y = -log10(`Pr(>|t|)`))

#Create the plot
l <- estlambda(res_out$`Pr(>|t|)`)
nom <- paste("CAD score vs. lipoprotein metabolites", "\n", "number of tests = ", nrow(res_out), "\nlambda = ", round(l$estimate, 3), ", se = ", round(l$se, 3), sep = "")
  p <- ggplot(res_out, aes(x = x, y = y, group = group, colour = subset)) +
    geom_point(aes(shape = factor(group)), size = 3) +
    geom_abline(slope = 1, intercept = 0, colour = "red") +
    #ggtitle(nom) +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 30), legend.text = element_text(size = 20)) +
    labs(x = expression(Expected ~ ~-log[10](P)), y = expression(Observed ~ ~-log[10](P))) #+
    #theme(legend.position = "right") +
    #scale_colour_manual(name = element_blank(),
    ###labels need to be in alphabetical order!!
    #labels = c("l-HDL", "l-VLDL", "LDL", "Rem P", "s-HDL", "vl-HDL"),
    #values = colours))

pdf("outputs/CAD_score_vs_lipoprotein_subclasses_QQ.pdf", width = 15, height = 10)
print(p)
dev.off()

#################################
# Differences between subsets ###
#################################
#Boxplot looking at effect size differences between groups
pdf("outputs/lipoprotein_subclass_effect_comparison.pdf", width = 15, height = 10)
ggplot(res_out) +
  geom_boxplot(aes(x = reorder(subset, abs(Estimate), FUN = median), y = abs(Estimate), fill = subset)) +
  labs(x = "Particle size class", y = "Relative effect estimate") +
  theme(text = element_text(size = 30), legend.position = "none", axis.title.x = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_discrete(labels=c("Small_HDL" = "Small HDL", "V_Large_HDL" = "Very large HDL", "Large_VLDL" = "Large VLDL", "Large_HDL" = "Large HDL", "Remnant_particles" = "Remnant particles", "LDL" = "LDL"))
dev.off()
#Testing for normality within group estimates
shap_test <- vector(mode = "list", length = length(unique(res_out$subset)))
names(shap_test) <- unique(res_out$subset)
for (i in res_out$subset) {
  temp_dat <- filter(res_out, subset == i)
  shap_test[[i]] <- shapiro.test(temp_dat$Estimate)
}
shap_test
#The estimates within many subsets are not normally distributed - use kruskal-wallis test

kw_test <- kruskal.test(Estimate ~ subset, data = res_out) 
dunn_test <- dunnTest(Estimate ~ subset, data = res_out, method = "bh")
dunn_res <- dunn_test[[2]]
write.table(dunn_res, "outputs/tables/lipoprotein_subclass_kw_test.txt", quote = F, col.names = T, row.names = F, sep = "\t")
################
# Associations #
################

CAD_assoc <- CAD_score_lr_nr %>%
  mutate(Bonferroni = p.adjust(`Pr(>|t|)`, method = "bonferroni")) %>%
  mutate(FDR = p.adjust(`Pr(>|t|)`, method = "fdr")) %>%
  filter(FDR < 0.05)

length(lipoproteins)  
sum(lipoproteins %in% CAD_assoc$Metabolite)

CAD_score_lr_nr <- arrange(CAD_score_lr_nr, `Pr(>|t|)`)
write.table(CAD_score_lr_nr, file = "outputs/tables/CAD_GRS_metabolite_full_assoc.txt", quote = F, col.names = T, row.names = F, sep = "\t")

#################
# subsets table #
#################
sub_tab <- res_out %>%
  arrange(subset) %>%
  select(Metabolite, subset, group)

write.table(sub_tab, file = "outputs/tables/subset_table.txt", quote = F, col.names = T, row.names = F, sep = "\t")

median(abs(CAD_score_lr_nr$Estimate))
print("CAD-GRS analysis complete - please proceed to the Individual variant analysis")





