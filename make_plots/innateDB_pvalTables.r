setwd(
  "~/Dropbox (ASU)/Placenta/BATCH2_PLACENTA_DECIDUA_ANALYSIS/HISAT_FeatureCounts/batch1_and_batch2/"
)
library(reshape)
library(ggpubr)
library(ggpubr)
library(dplyr)
library(lattice)
library(car)
library(splitstackshape)

df_merged_pisarska <-
  read.delim ("genelists/pisarska/cpm_placenta_df_genes.txt")
df_merged_pisarska$group <- "late first trimester"
df_merged_placenta <-
  read.delim ("genelists/placentas_batch_1and2/cpm_placenta_df_genes.txt")
df_merged_placenta$group <- "term >= 36 weeks"
df_decidua <-read.delim ("genelists/deciduas/cpm_decidua_df_genes.txt")
df_decidua$group <- "decidua"

setwd(
  "~/Dropbox (ASU)/Placenta/BATCH2_PLACENTA_DECIDUA_ANALYSIS/HISAT_FeatureCounts/batch1_and_batch2/innateDB/"
)

innateDB <- read.delim("innateDB.txt", header = TRUE)
genes <- read.csv("~/Dropbox (ASU)/Placenta/BATCH2_PLACENTA_DECIDUA_ANALYSIS/HISAT_FeatureCounts/batch1_and_batch2/counts_pheno/genesID.csv", header = TRUE)
innateDB_Human <- intersect(innateDB$Geneid, genes$Geneid)
rbind.all.columns <- function(x, y) {
  x.diff <- setdiff(colnames(x), colnames(y))
  y.diff <- setdiff(colnames(y), colnames(x))
  x[, c(as.character(y.diff))] <- NA
  y[, c(as.character(x.diff))] <- NA
  return(rbind(x, y))
}
df_merged <- rbind.all.columns(df_merged_placenta, df_merged_pisarska)
df_merged$tissue <- c("placenta")
df_merged_placentaAndDecidua <- rbind.all.columns(df_merged, df_decidua)

# gene lists 
innateDB_inter_pisarska <- intersect(innateDB$Geneid, df_merged_pisarska$Geneid)
innateDB_inter_placenta <- intersect(innateDB$Geneid, df_merged_placenta$Geneid)
innateDB_inter_decidua <- intersect(innateDB$Geneid, df_decidua$Geneid)
FirstAndTermGenes <- intersect(innateDB_inter_pisarska, innateDB_inter_placenta)
placentaAndDeciduaGenes <- intersect(innateDB_inter_decidua, FirstAndTermGenes)
#----------------------------------- 
# If the variances of the two groups being compared are different (heteroscedasticity), 
# it’s possible to use the Welch t test, an adaptation of Student t-test.
# t.test(x, y, alternative = "two.sided", var.equal = FALSE)
# x,y: numeric vectors
# alternative: the alternative hypothesis. Allowed value is one of “two.sided” (default), “greater” or “less”.
# var.equal: a logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch test is used.
#  test for homogeneity in variances.

# female vs male expression in First trimester 
df_pvals_lateFirst <- data.frame()
vt_pval <- NULL
ttest_pval <- NULL
sha_f_pval <- NULL
sha_m_pval <- NULL
wil_pval <- NULL
female_mean <- NULL
male_mean <- NULL
female_median <- NULL
male_median <- NULL
Geneid <- NULL
for(i in innateDB_inter_pisarska){
  Geneid <- c(Geneid, i)
  geneDF <- subset(df_merged_pisarska, Geneid == i)
  # geneDF <- subset(geneDF, group == "late first trimester")
  # are the variance equal between the groups
  vt <- var.test(value ~sex, geneDF)
  vt_pval <- c(vt_pval, vt$p.value)
  # are the means equal between the two unpaired groups
  ttest <- t.test(value ~sex, geneDF)#, var.equal = TRUE)
  ttest_pval <- c(ttest_pval, ttest$p.value)
  # is the data in each group normally distributed
  # if the data is NOT normally distributed, 
  # use a Wilcox rank sum test
  wil <- wilcox.test(value ~sex, geneDF)
  wil_pval <- c(wil_pval, wil$p.value)
  df_F <- subset(geneDF, sex == "female")
  female_mean <- c(female_mean, mean(df_F$value))
  female_median <- c(female_median, median(df_F$value))
  sha_f <- shapiro.test(df_F$value)
  sha_f_pval <- c(sha_f_pval, sha_f$p.value)
  df_M <- subset(geneDF, sex == "male")
  male_mean <- c(male_mean, mean(df_M$value))
  male_median <- c(male_median, median(df_M$value))
  sha_m <- shapiro.test(df_M$value)
  sha_m_pval <- c(sha_m_pval, sha_m$p.value)
  df_pvals_lateFirst <- cbind(Geneid, female_mean, female_median, male_mean, male_median, vt_pval, sha_f_pval, sha_m_pval, ttest_pval, wil_pval)
}
head(df_pvals_lateFirst)
write.table(df_pvals_lateFirst, "df_pvals_lateFirst_ImmuneDBgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

df_pvals_lateFirst <- read.delim("df_pvals_lateFirst_ImmuneDBgenes.txt", header = TRUE, sep = "\t")
df_pvals_lateFirst_ImmuneDBgenes_sig05 <- subset(df_pvals_lateFirst, wil_pval <= 0.05)
write.table(df_pvals_lateFirst_ImmuneDBgenes_sig05, "df_pvals_lateFirst_ImmuneDBgenes_sig05.txt", sep = "\t", quote = FALSE, row.names = FALSE)

df_pvals_lateFirst_ImmuneDBgenes_sig05_femaleBias <- subset(df_pvals_lateFirst_ImmuneDBgenes_sig05, female_median > male_median) 
df_pvals_lateFirst_ImmuneDBgenes_sig05_maleBias <- subset(df_pvals_lateFirst_ImmuneDBgenes_sig05, female_median < male_median) 
write.table(df_pvals_lateFirst_ImmuneDBgenes_sig05_femaleBias, "df_pvals_lateFirst_ImmuneDBgenes_sig05_femaleBias.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(df_pvals_lateFirst_ImmuneDBgenes_sig05_maleBias, "df_pvals_lateFirst_ImmuneDBgenes_sig05_maleBias.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# female vs male expression in term Wilson placentas
df_pvals_term <- data.frame()
vt_pval <- NULL
ttest_pval <- NULL
sha_f_pval <- NULL
sha_m_pval <- NULL
wil_pval <- NULL
female_mean <- NULL
male_mean <- NULL
female_median <- NULL
male_median <- NULL
Geneid <- NULL
for(i in innateDB_inter_placenta){
  Geneid <- c(Geneid, i)
  geneDF <- subset(df_merged_placenta, Geneid == i)
  # geneDF <- subset(geneDF, group == "late first trimester")
  # are the variance equal between the groups
  vt <- var.test(value ~sex, geneDF)
  vt_pval <- c(vt_pval, vt$p.value)
  # are the means equal between the two unpaired groups
  ttest <- t.test(value ~sex, geneDF)#, var.equal = TRUE)
  ttest_pval <- c(ttest_pval, ttest$p.value)
  # is the data in each group normally distributed
  # if the data is NOT normally distributed, 
  # use a Wilcox rank sum test
  wil <- wilcox.test(value ~sex, geneDF)
  wil_pval <- c(wil_pval, wil$p.value)
  df_F <- subset(geneDF, sex == "female")
  female_mean <- c(female_mean, mean(df_F$value))
  female_median <- c(female_median, median(df_F$value))
 # sha_f <- shapiro.test(df_F$value)
 # sha_f_pval <- c(sha_f_pval, sha_f$p.value)
  df_M <- subset(geneDF, sex == "male")
  male_mean <- c(male_mean, mean(df_M$value))
  male_median <- c(male_median, median(df_M$value))
  # sha_m <- shapiro.test(as.numeric(df_M$value))
  # sha_m_pval <- c(sha_m_pval, sha_m$p.value)
  df_pvals_term <- cbind(Geneid, female_mean, female_median, male_mean, male_median, vt_pval, ttest_pval, wil_pval)
}
df_pvals_term
write.table(df_pvals_term, "df_pvals_term_ImmuneDBgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

df_pvals_term_ImmuneDBgenes_all <- read.delim("df_pvals_term_ImmuneDBgenes.txt", header = TRUE, sep = "\t")
df_pvals_term_ImmuneDBgenes_sig05 <- subset(df_pvals_term_ImmuneDBgenes_all, wil_pval <= 0.05)
write.table(df_pvals_term_ImmuneDBgenes_sig05, "df_pvals_term_ImmuneDBgenes_sig05.txt", sep = "\t", quote = FALSE, row.names = FALSE)

df_pvals_term_ImmuneDBgenes_sig05_femaleBias <- subset(df_pvals_term_ImmuneDBgenes_sig05, female_median > male_median) 
df_pvals_term_ImmuneDBgenes_sig05_maleBias <- subset(df_pvals_term_ImmuneDBgenes_sig05, female_median < male_median) 
write.table(df_pvals_term_ImmuneDBgenes_sig05_femaleBias, "df_pvals_term_ImmuneDBgenes_sig05_femaleBias.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(df_pvals_term_ImmuneDBgenes_sig05_maleBias, "df_pvals_term_ImmuneDBgenes_sig05_maleBias.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# # female vs male expression in ALL placentas, both first and term 
# df_pvals_FemaleVsMale_allPlacentas <- data.frame()
# vt_pval <- NULL
# ttest_pval <- NULL
# sha_f_pval <- NULL
# sha_m_pval <- NULL
# wil_pval <- NULL
# female_mean <- NULL
# male_mean <- NULL
# female_median <- NULL
# male_median <- NULL
# Geneid <- NULL
# for(i in innateDB_inter_pisarska){
#   Geneid <- c(Geneid, i)
#   geneDF <- subset(df_merged, Geneid == i)
#   # geneDF <- subset(geneDF, group == "late first trimester")
#   # are the variance equal between the groups
#   vt <- var.test(value ~sex, geneDF)
#   vt_pval <- c(vt_pval, vt$p.value)
#   # are the means equal between the two unpaired groups
#   ttest <- t.test(value ~sex, geneDF)#, var.equal = TRUE)
#   ttest_pval <- c(ttest_pval, ttest$p.value)
#   # is the data in each group normally distributed
#   # if the data is NOT normally distributed, 
#   # use a Wilcox rank sum test
#   wil <- wilcox.test(value ~sex, geneDF)
#   wil_pval <- c(wil_pval, wil$p.value)
#   df_F <- subset(geneDF, sex == "female")
#   female_mean <- c(female_mean, mean(df_F$value))
#   female_median <- c(female_median, median(df_F$value))
#   sha_f <- shapiro.test(df_F$value)
#   sha_f_pval <- c(sha_f_pval, sha_f$p.value)
#   df_M <- subset(geneDF, sex == "male")
#   male_mean <- c(male_mean, mean(df_M$value))
#   male_median <- c(male_median, median(df_M$value))
#   sha_m <- shapiro.test(df_M$value)
#   sha_m_pval <- c(sha_m_pval, sha_m$p.value)
#   df_pvals_FemaleVsMale_allPlacentas <- cbind(Geneid, female_mean, female_median, male_mean, male_median, vt_pval, sha_f_pval, sha_m_pval, ttest_pval, wil_pval)
# }
# write.table(df_pvals_FemaleVsMale_allPlacentas, "df_pvals_FemaleVsMale_allPlacentas_ImmuneDBgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# 
# df_pvals_FemaleVsMale_allPlacentas_ImmuneDBgenes_all <- read.delim("df_pvals_FemaleVsMale_allPlacentas_ImmuneDBgenes.txt", header = TRUE, sep = "\t")
# df_pvals_FemaleVsMale_allPlacentas_ImmuneDBgenes_sig05 <- subset(df_pvals_FemaleVsMale_allPlacentas_ImmuneDBgenes_all, wil_pval <= 0.05)
# write.table(df_pvals_FemaleVsMale_allPlacentas_ImmuneDBgenes_sig05, "df_pvals_FemaleVsMale_allPlacentas_ImmuneDBgenes_sig05.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# 
# df_pvals_FemaleVsMale_allPlacentas_ImmuneDBgenes_sig05_femaleBias <- subset(df_pvals_FemaleVsMale_allPlacentas_ImmuneDBgenes_sig05, female_median > male_median) 
# df_pvals_FemaleVsMale_allPlacentas_ImmuneDBgenes_sig05_maleBias <- subset(df_pvals_FemaleVsMale_allPlacentas_ImmuneDBgenes_sig05, female_median < male_median) 
# write.table(df_pvals_FemaleVsMale_allPlacentas_ImmuneDBgenes_sig05_femaleBias, "df_pvals_FemaleVsMale_allPlacentas_ImmuneDBgenes_sig05_femaleBias.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# write.table(df_pvals_FemaleVsMale_allPlacentas_ImmuneDBgenes_sig05_maleBias, "df_pvals_FemaleVsMale_allPlacentas_ImmuneDBgenes_sig05_maleBias.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# first vs term expression 
df_pvals_firstVsterm_allPlacentas <- data.frame()
vt_pval <- NULL
ttest_pval <- NULL
sha_f_pval <- NULL
sha_m_pval <- NULL
wil_pval <- NULL
first_mean <- NULL
term_mean <- NULL
first_median <- NULL
term_median <- NULL
Geneid <- NULL

for(i in FirstAndTermGenes){
  Geneid <- c(Geneid, i)
  geneDF <- subset(df_merged, Geneid == i)
  # geneDF <- subset(geneDF, group == "late first trimester")
  # are the variance equal between the groups
  vt <- var.test(value ~group, geneDF)
  vt_pval <- c(vt_pval, vt$p.value)
  # are the means equal between the two unpaired groups
  ttest <- t.test(value ~group, geneDF)#, var.equal = TRUE)
  ttest_pval <- c(ttest_pval, ttest$p.value)
  # is the data in each group normally distributed
  # if the data is NOT normally distributed, 
  # use a Wilcox rank sum test
  wil <- wilcox.test(value ~sex, geneDF)
  wil_pval <- c(wil_pval, wil$p.value)
  df_F <- subset(geneDF, group == "late first trimester")
  first_mean <- c(first_mean, mean(df_F$value))
  first_median <- c(first_median, median(df_F$value))
  sha_f <- shapiro.test(df_F$value)
  sha_f_pval <- c(sha_f_pval, sha_f$p.value)
  df_M <- subset(geneDF, group == "term >= 36 weeks")
  term_mean <- c(term_mean, mean(df_M$value))
  term_median <- c(term_median, median(df_M$value))
  sha_m <- shapiro.test(df_M$value)
  sha_m_pval <- c(sha_m_pval, sha_m$p.value)
  df_pvals_firstVsterm_allPlacentas <- cbind(Geneid, first_mean, first_median, term_mean, term_median, vt_pval, sha_f_pval, sha_m_pval, ttest_pval, wil_pval)
}
write.table(df_pvals_firstVsterm_allPlacentas, "df_pvals_firstVsterm_allPlacentas_ImmuneDBgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

df_pvals_firstVsterm_allPlacentas_ImmuneDBgenes_all <- read.delim("df_pvals_firstVsterm_allPlacentas_ImmuneDBgenes.txt", header = TRUE, sep = "\t")
df_pvals_firstVsterm_allPlacentas_ImmuneDBgenes_sig05 <- subset(df_pvals_firstVsterm_allPlacentas_ImmuneDBgenes_all, wil_pval <= 0.05)
write.table(df_pvals_firstVsterm_allPlacentas_ImmuneDBgenes_sig05, "df_pvals_firstVsterm_allPlacentas_ImmuneDBgenes_sig05.txt", sep = "\t", quote = FALSE, row.names = FALSE)

df_pvals_firstVsterm_allPlacentas_ImmuneDBgenes_sig05_firstBias <- subset(df_pvals_firstVsterm_allPlacentas_ImmuneDBgenes_sig05, first_median > term_median) 
df_pvals_firstVsterm_allPlacentas_ImmuneDBgenes_sig05_termBias <- subset(df_pvals_firstVsterm_allPlacentas_ImmuneDBgenes_sig05, first_median < term_median) 
write.table(df_pvals_firstVsterm_allPlacentas_ImmuneDBgenes_sig05_firstBias, "df_pvals_firstVsterm_allPlacentas_ImmuneDBgenes_sig05_firstBias.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(df_pvals_firstVsterm_allPlacentas_ImmuneDBgenes_sig05_termBias, "df_pvals_firstVsterm_allPlacentas_ImmuneDBgenes_sig05_termBias.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# decidua - male vs female offspring sex 
df_pvals_decidua <- data.frame()
vt_pval <- NULL
ttest_pval <- NULL
sha_f_pval <- NULL
sha_m_pval <- NULL
wil_pval <- NULL
female_mean <- NULL
male_mean <- NULL
female_median <- NULL
male_median <- NULL
Geneid <- NULL
for(i in innateDB_inter_decidua){
  Geneid <- c(Geneid, i)
  geneDF <- subset(df_decidua, Geneid == i)
  # are the variance equal between the groups
  vt <- var.test(value ~Offspring_sex, geneDF)
  vt_pval <- c(vt_pval, vt$p.value)
  # are the means equal between the two unpaired groups
  ttest <- t.test(value ~Offspring_sex, geneDF)#, var.equal = TRUE)
  ttest_pval <- c(ttest_pval, ttest$p.value)
  # is the data in each group normally distributed
  # if the data is NOT normally distributed, 
  # use a Wilcox rank sum test
  wil <- wilcox.test(value ~Offspring_sex, geneDF)
  wil_pval <- c(wil_pval, wil$p.value)
  df_F <- subset(geneDF, Offspring_sex == "female")
  female_mean <- c(female_mean, mean(df_F$value))
  female_median <- c(female_median, median(df_F$value))
  sha_f <- shapiro.test(df_F$value)
  sha_f_pval <- c(sha_f_pval, sha_f$p.value)
  df_M <- subset(geneDF, Offspring_sex == "male")
  male_mean <- c(male_mean, mean(df_M$value))
  male_median <- c(male_median, median(df_M$value))
  sha_m <- shapiro.test(df_M$value)
  sha_m_pval <- c(sha_m_pval, sha_m$p.value)
  df_pvals_decidua <- cbind(Geneid, female_mean, female_median, male_mean, male_median, vt_pval, sha_f_pval, sha_m_pval, ttest_pval, wil_pval)
}
write.table(df_pvals_decidua, "df_pvals_decidua_ImmuneDBgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

df_pvals_decidua_ImmuneDBgenes_all <- read.delim("df_pvals_decidua_ImmuneDBgenes.txt", header = TRUE, sep = "\t")
df_pvals_decidua_ImmuneDBgenes_sig05 <- subset(df_pvals_decidua_ImmuneDBgenes_all, wil_pval <= 0.05)
write.table(df_pvals_decidua_ImmuneDBgenes_sig05, "df_pvals_decidua_ImmuneDBgenes_sig05.txt", sep = "\t", quote = FALSE, row.names = FALSE)

df_pvals_decidua_ImmuneDBgenes_sig05_femaleBias <- subset(df_pvals_decidua_ImmuneDBgenes_sig05, female_median > male_median) 
df_pvals_decidua_ImmuneDBgenes_sig05_maleBias <- subset(df_pvals_decidua_ImmuneDBgenes_sig05, female_median < male_median) 
write.table(df_pvals_decidua_ImmuneDBgenes_sig05_femaleBias, "df_pvals_decidua_ImmuneDBgenes_sig05_femaleBias.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(df_pvals_decidua_ImmuneDBgenes_sig05_maleBias, "df_pvals_decidua_ImmuneDBgenes_sig05_maleBias.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# placenta vs decidua expression 
df_pvals_placentaVsdecidua_allPlacentas <- data.frame()
vt_pval <- NULL
ttest_pval <- NULL
sha_f_pval <- NULL
sha_m_pval <- NULL
wil_pval <- NULL
placenta_mean <- NULL
decidua_mean <- NULL
placenta_median <- NULL
decidua_median <- NULL
Geneid <- NULL
for(i in placentaAndDeciduaGenes){
  Geneid <- c(Geneid, i)
  geneDF <- subset(df_merged_placentaAndDecidua, Geneid == i)
  # are the variance equal between the groups
  vt <- var.test(value ~tissue, geneDF)
  vt_pval <- c(vt_pval, vt$p.value)
  # are the means equal between the two unpaired groups
  ttest <- t.test(value ~tissue, geneDF)#, var.equal = TRUE)
  ttest_pval <- c(ttest_pval, ttest$p.value)
  # is the data in each group normally distributed
  # if the data is NOT normally distributed, 
  # use a Wilcox rank sum test
  wil <- wilcox.test(value ~tissue, geneDF)
  wil_pval <- c(wil_pval, wil$p.value)
  df_F <- subset(geneDF, tissue == "placenta")
  placenta_mean <- c(placenta_mean, mean(df_F$value))
  placenta_median <- c(placenta_median, median(df_F$value))
  sha_f <- shapiro.test(df_F$value)
  sha_f_pval <- c(sha_f_pval, sha_f$p.value)
  df_M <- subset(geneDF, tissue == "decidua")
  decidua_mean <- c(decidua_mean, mean(df_M$value))
  decidua_median <- c(decidua_median, median(df_M$value))
  sha_m <- shapiro.test(df_M$value)
  sha_m_pval <- c(sha_m_pval, sha_m$p.value)
  df_pvals_placentaVsdecidua_allPlacentas <- cbind(Geneid, placenta_mean, placenta_median, decidua_mean, decidua_median, vt_pval, sha_f_pval, sha_m_pval, ttest_pval, wil_pval)
}
write.table(df_pvals_placentaVsdecidua_allPlacentas, "df_pvals_placentaVsdecidua_allPlacentas_ImmuneDBgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

df_pvals_placentaVsdecidua_allPlacentas_ImmuneDBgenes_all <- read.delim("df_pvals_placentaVsdecidua_allPlacentas_ImmuneDBgenes.txt", header = TRUE, sep = "\t")
df_pvals_placentaVsdecidua_allPlacentas_ImmuneDBgenes_sig05 <- subset(df_pvals_placentaVsdecidua_allPlacentas_ImmuneDBgenes_all, wil_pval <= 0.05)
write.table(df_pvals_placentaVsdecidua_allPlacentas_ImmuneDBgenes_sig05, "df_pvals_placentaVsdecidua_allPlacentas_ImmuneDBgenes_sig05.txt", sep = "\t", quote = FALSE, row.names = FALSE)

df_pvals_placentaVsdecidua_allPlacentas_ImmuneDBgenes_sig05_placentaBias <- subset(df_pvals_placentaVsdecidua_allPlacentas_ImmuneDBgenes_sig05, placenta_median > decidua_median) 
df_pvals_placentaVsdecidua_allPlacentas_ImmuneDBgenes_sig05_deciduaBias <- subset(df_pvals_placentaVsdecidua_allPlacentas_ImmuneDBgenes_sig05, placenta_median < decidua_median) 
write.table(df_pvals_placentaVsdecidua_allPlacentas_ImmuneDBgenes_sig05_placentaBias, "df_pvals_placentaVsdecidua_allPlacentas_ImmuneDBgenes_sig05_placentaBias.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(df_pvals_placentaVsdecidua_allPlacentas_ImmuneDBgenes_sig05_deciduaBias, "df_pvals_placentaVsdecidua_allPlacentas_ImmuneDBgenes_sig05_deciduaBias.txt", sep = "\t", quote = FALSE, row.names = FALSE)


