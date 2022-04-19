setwd("/Users/m239830/Dropbox (ASU)/Placenta/BATCH2_PLACENTA_DECIDUA_ANALYSIS/HISAT_FeatureCounts/batch1_and_batch2/GTEx_TPM/GTEx_MeanTPM")
library(stringr)
library(tidyverse)
library(dplyr)
#install.packages("corrplot")
source("http://www.sthda.com/upload/rquery_cormat.r")
library(corrplot)
library(UpSetR)
library(ggplot2)
library(dplyr)
library(Hmisc)  

GTEx <- c("Skin-SunExposed_Lowerleg",
           "WholeBlood",
           "Thyroid",
           "Adipose-Subcutaneous",
           "Artery-Tibial",
           "Skin-NotSunExposed_Suprapubic",
           "Lung",
           "Esophagus-Muscularis",
           "Colon-Transverse",
           "Breast-MammaryTissue",
           "Adipose-Visceral_Omentum",
           "Heart-LeftVentricle",
           "Artery-Aorta",
           "Brain-Cortex",
           "Heart-AtrialAppendage",
           "Cells-Culturedfibroblasts",
           "Esophagus-GastroesophagealJunction",
           "Stomach",
           "Pancreas",
           "Pituitary",
           "Brain-Caudate_basalganglia",
           "AdrenalGland",
           "Esophagus-Mucosa",
           "Artery-Coronary",
           "Spleen",
           "Brain-Hippocampus",
           "Brain-Cerebellum",
           "Brain-Hypothalamus",
           "Brain-CerebellarHemisphere",
           "Brain-Putamen_basalganglia",
           "Brain-FrontalCortex_BA9",
           "Colon-Sigmoid",
           "MinorSalivaryGland",
           "Cells-EBV-transformedlymphocytes",
           "Liver",
           "Brain-Spinalcord_cervicalc-1",
           "Brain-Amygdala",
           "Brain-substantianigra",
           "Muscle-Skeletal",
           "SmallIntestine-TerminalIleum",
           "Kidney-Cortex",
           "Brain-Nucleusaccumbens_basalganglia",
           "Nerve-Tibial",
           "Brain-Anteriorcingulatecortex_BA24",
           "Bladder")
# read in  placenta tissue files
term_placenta <- read.delim("/Users/m239830/Dropbox (ASU)/Placenta/BATCH2_PLACENTA_DECIDUA_ANALYSIS/HISAT_FeatureCounts/batch1_and_batch2/DEGs/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/DEGs_placentas_genesFvsM_fdr1_lfc0.txt", sep = '\t', header = TRUE)
firstTrimester_placenta <- read.delim("/Users/m239830/Dropbox (ASU)/Placenta/BATCH2_PLACENTA_DECIDUA_ANALYSIS/HISAT_FeatureCounts/batch1_and_batch2/DEGs/pisarska/DEGs_placentas_genesFvsM_fdr1_lfc0.txt", sep = '\t', header = TRUE)
placentaSexGenes <- subset(term_placenta, chr == "chrX" | chr =="chrY")

for(x in GTEx) {
  t <- read.delim(paste0(x, "_bysex_counts_mean.txt")) # load file
  t$logFC <- log2((t$female_mean+.01)/(t$male_mean+.01))
  names(t)[names(t) == "Description"] <- "Geneid"
  tissue <- subset(t, (t$Geneid %in% term_placenta$Geneid))
  assign(paste(x, sep = ''), tissue)
}
#subset(firstTrimester_placenta, firstTrimester_placenta$Geneid == "MANEA")
All_Tissues <- cbind(`Skin-SunExposed_Lowerleg`,
                     `WholeBlood`,
                     `Thyroid`,
                     `Adipose-Subcutaneous`,
                     `Artery-Tibial`,
                     `Skin-NotSunExposed_Suprapubic`,
                     `Lung`,
                     `Esophagus-Muscularis`,
                     `Colon-Transverse`,
                     `Breast-MammaryTissue`,
                     `Adipose-Visceral_Omentum`,
                     `Heart-LeftVentricle`,
                     `Artery-Aorta`,
                     `Brain-Cortex`,
                     `Heart-AtrialAppendage`,
                     #`Cells-Culturedfibroblasts`,
                     `Esophagus-GastroesophagealJunction`,
                     `Stomach`,
                     `Pancreas`,
                     `Pituitary`,
                     `Brain-Caudate_basalganglia`,
                     `AdrenalGland`,
                     `Esophagus-Mucosa`,
                     `Artery-Coronary`,
                     `Spleen`,
                     `Brain-Hippocampus`,
                     `Brain-Cerebellum`,
                     `Brain-Hypothalamus`,
                     `Brain-CerebellarHemisphere`,
                     `Brain-Putamen_basalganglia`,
                     `Brain-FrontalCortex_BA9`,
                     `Colon-Sigmoid`,
                     `MinorSalivaryGland`,
                     #`Cells-EBV-transformedlymphocytes`,
                     `Liver`,
                     `Brain-Spinalcord_cervicalc-1`,
                     `Brain-Amygdala`,
                     `Brain-substantianigra`,
                   #  `Muscle-Skeletal`,
                     `SmallIntestine-TerminalIleum`,
                     `Kidney-Cortex`,
                     `Brain-Nucleusaccumbens_basalganglia`,
                     `Nerve-Tibial`,
                     `Brain-Anteriorcingulatecortex_BA24`,
                     `Bladder`)
drops <- c("Name", "Geneid", "female_mean","male_mean")
GTEx_log2FC <- All_Tissues[ ,!(names(All_Tissues) %in% drops)]
GTEx_log2FC$Geneid <- WholeBlood$Geneid


term_placenta_df <- data.frame(term_placenta$Geneid, term_placenta$logFC)
firstTrimester_placenta_df <- data.frame(firstTrimester_placenta$Geneid, firstTrimester_placenta$logFC)
colnames(term_placenta_df) <- c("Geneid", "logFC")
colnames(firstTrimester_placenta_df) <- c("Geneid", "logFC")

MergeTissuesByID <- list(`term_placenta_df`,
     `firstTrimester_placenta_df`,
     `GTEx_log2FC`) %>% reduce(full_join, by = "Geneid")


new <- MergeTissuesByID[!duplicated(MergeTissuesByID[2]),]
rownames(new) <- new[, 1]
df_tissues <-
  new[,grep("logFC", colnames(new))]

colnames(df_tissues) <- c("Term placenta",
                          "Late first trimester placenta",
                          "Skin sun exposed lower leg",
                          "Whole blood",
                          "Thyroid",
                          "Adipose subcutaneous",
                          "Artery tibial",
                          "Skin not sun exposed suprapubic",
                          "Lung",
                          "Esophagus muscularis",
                          "Colon transverse",
                          "Breast mammary",
                          "Adipose visceral omentum",
                          "Heart left ventricle",
                          "Artery aorta",
                          "Brain cortex",
                          "Heart atrial appendage",
                        #  "Cells-Culturedfibroblasts",
                          "Esophagus",
                          "Stomach",
                          "Pancreas",
                          "Pituitary",
                          "Brain caudate",
                          "Adrenal gland",
                          "Esophagus mucosa",
                          "Artery coronary",
                          "Spleen",
                          "Brain hippocampus",
                          "Brain cerebellum",
                          "Brain hypothalamus",
                          "Brain cerebellar hemisphere",
                          "Brain putamen",
                          "Brain frontal cortex",
                          "Colon",
                          "Minor salivary gland",
                        #  "Cells-EBV-transformedlymphocytes",
                          "Liver",
                          "Brain spinal cord",
                          "Brain amygdala",
                          "Brain substantia nigra",
                          #  "Muscle-Skeletal",
                          "Small intestine",
                          "Kidney",
                          "Brain nucleus accumbens",
                          "Nerve tibial",
                          "Brain anterior cingulate cortex",
                          "Bladder")
df_tissues_noNA <- na.omit(df_tissues)
#---------------
# correlation - all genes 
#---------------
tissueMatrix <- as.matrix(df_tissues_noNA)
# dim(tissueMatrix)
# sexDiffGenes <- subset(term_placenta, term_placenta$adj.P.Val < 0.05)
# tissueMatrix <- subset(tissueMatrix, !(row.names(tissueMatrix) %in% sexDiffGenes$Geneid))
M <- cor(tissueMatrix)
M_df <- as.data.frame(M)
M_df$`Term placenta` <- as.numeric(M_df$`Term placenta`)
Mdf_round <- M_df %>% mutate_if(is.numeric, round, digits=2)
Mdf_round$Tissues <- row.names(M_df)
Mdf_round$Tissues <- factor(Mdf_round$Tissues, levels = Mdf_round$Tissues[order(Mdf_round$`Term placenta`)])
# plot 
p <- ggplot(Mdf_round, aes(x = 1,y = Tissues, fill = `Term placenta`)) + 
  scale_fill_gradient(low = "white", high = "steelblue") +
  geom_tile() + 
  ylab("")

P2<- p +  geom_text(aes(label=`Term placenta`)) +  theme(axis.title.x=element_blank(),
                                                  axis.text.x=element_blank(),
                                                  axis.ticks.x=element_blank())

P2
ggsave("r2_allGenes.pdf", 
       P2,
       width = 5, height = 7)
dev.off()
dev.off()

#-----
# all genes, sex chromosome linked genes, autosomal genes
#-----
tissueMatrix <- as.matrix(df_tissues_noNA)
# dim(tissueMatrix)
# sexDiffGenes <- subset(term_placenta, term_placenta$adj.P.Val < 0.05)
# tissueMatrix <- subset(tissueMatrix, !(row.names(tissueMatrix) %in% sexDiffGenes$Geneid))
M <- cor(tissueMatrix)
M_test <- rcorr(tissueMatrix)

#M_sexChr_test$P
M_testP <- as.data.frame(M_test$P)
M_testP <- as.numeric(M_testP$`Term placenta`)
M_df <- as.data.frame(M)
M_df$`Term placenta` <- as.numeric(M_df$`Term placenta`)
Mdf_round <- M_df %>% mutate_if(is.numeric, round, digits=3)
Mdf_round$Tissues <- row.names(M_df)
Mdf_round$Tissues <- factor(Mdf_round$Tissues, levels = Mdf_round$Tissues[order(Mdf_round$`Term placenta`)])

dim(tissueMatrix)
tissueMatrix <- as.matrix(df_tissues_noNA)

# chr X
sexChrGenes_X <- subset(term_placenta, chr == "chrX")
sexChrGenes_X_tissueMatrix <- subset(tissueMatrix, (row.names(tissueMatrix) %in% sexChrGenes_X$Geneid))
dim(sexChrGenes_X_tissueMatrix)
M_sexChr_X <- cor(sexChrGenes_X_tissueMatrix)
M_sexChr_X_test <- rcorr(sexChrGenes_X_tissueMatrix)
#M_sexChr_test$P
M_sexChr_X_testP <- as.data.frame(M_sexChr_X_test$P)
M_sexChr_X_testP <- as.numeric(M_sexChr_X_testP$`Term placenta`)
M_sexChr_X_df <- as.data.frame(M_sexChr_X)
M_sexChr_X_df$`Term placenta` <- as.numeric(M_sexChr_X_df$`Term placenta`)
M_sexChr_X_df_round <- M_sexChr_X_df %>% mutate_if(is.numeric, round, digits=3)
M_sexChr_X_df_round$Tissues <- row.names(M_sexChr_X_df)
M_sexChr_X_df_round$Tissues <- factor(M_sexChr_X_df_round$Tissues, levels = M_sexChr_X_df_round$Tissues[order(M_sexChr_X_df_round$`Term placenta`)])

# chr Y
sexChrGenes_Y <- subset(term_placenta, chr == "chrY")
sexChrGenes_Y_tissueMatrix <- subset(tissueMatrix, (row.names(tissueMatrix) %in% sexChrGenes_Y$Geneid))
dim(sexChrGenes_Y_tissueMatrix)
M_sexChr_Y <- cor(sexChrGenes_Y_tissueMatrix)
M_sexChr_Y_test <- rcorr(sexChrGenes_Y_tissueMatrix)
#M_sexChr_test$P
M_sexChr_Y_testP <- as.data.frame(M_sexChr_Y_test$P)
M_sexChr_Y_testP <- as.numeric(M_sexChr_Y_testP$`Term placenta`)
M_sexChr_Y_df <- as.data.frame(M_sexChr_Y)
M_sexChr_Y_df$`Term placenta` <- as.numeric(M_sexChr_Y_df$`Term placenta`)
M_sexChr_Y_df_round <- M_sexChr_Y_df %>% mutate_if(is.numeric, round, digits=3)
M_sexChr_Y_df_round$Tissues <- row.names(M_sexChr_Y_df)
M_sexChr_Y_df_round$Tissues <- factor(M_sexChr_Y_df_round$Tissues, levels = M_sexChr_Y_df_round$Tissues[order(M_sexChr_Y_df_round$`Term placenta`)])

#-------------
tissueMatrix <- as.matrix(df_tissues_noNA)
autosomal_all <- subset(term_placenta, chr != "chrX" & chr != "chrY")
autosomal_all_tissueMatrix <- subset(tissueMatrix, (row.names(tissueMatrix) %in% autosomal_all$Geneid))
dim(autosomal_all_tissueMatrix)
M_autosomal_all <- cor(autosomal_all_tissueMatrix)
M_autosomal_all_test <- rcorr(M_autosomal_all)
#M_autosomal_all_test$P
M_autosomal_all_testP <- as.data.frame(M_autosomal_all_test$P)
M_autosomal_all_testP_rownames <- rownames(M_autosomal_all_testP)
M_autosomal_all_testP <- as.numeric(M_autosomal_all_testP$`Term placenta`)
df_M_autorownames <- as.data.frame(M_autosomal_all_testP_rownames)
df_M_autosomal_all_testP <- as.data.frame(M_autosomal_all_testP)
rownames(df_M_autosomal_all_testP) <- df_M_autorownames$M_autosomal_all_testP_rownames
M_autosomal_all_df <- as.data.frame(M_autosomal_all)
M_autosomal_all_df$`Term placenta` <- as.numeric(M_autosomal_all_df$`Term placenta`)
M_autosomal_all_df_round <- M_autosomal_all_df %>% mutate_if(is.numeric, round, digits=3)
M_autosomal_all_df_round$Tissues <- row.names(M_autosomal_all_df)
M_autosomal_all_df_round$Tissues <- factor(M_autosomal_all_df_round$Tissues, levels = M_autosomal_all_df_round$Tissues[order(M_autosomal_all_df_round$`Term placenta`)])

#--- all genes, sex-chromosome linked and autosomal 
Mdf_round$corr <- "All genes"
#M_sexChr_df_round$corr <- "Sex-linked"
M_autosomal_all_df_round$corr <- "Autosomal"
M_sexChr_X_df_round$corr <- "X-linked"
M_sexChr_Y_df_round$corr <- "Y-linked"


groupCorr <- rbind(Mdf_round, M_sexChr_X_df_round, M_sexChr_Y_df_round, M_autosomal_all_df_round)
library(reshape2)
meltCorr <- melt(groupCorr)
termCorr <- subset(meltCorr, meltCorr$variable == "Term placenta")
termCorr$pvalue <- NULL

# Pvalues 
M_testP_df <- data.frame(M_testP)
names(M_testP_df)[names(M_testP_df) == "M_testP"] <- "pvalue"
# M_sexChr_testP_df <- data.frame(M_sexChr_testP)
# names(M_sexChr_testP_df)[names(M_sexChr_testP_df) == "M_sexChr_testP"] <- "pvalue"
M_autosomal_all_testP_df <- data.frame(M_autosomal_all_testP)
names(M_autosomal_all_testP_df)[names(M_autosomal_all_testP_df) == "M_autosomal_all_testP"] <- "pvalue"
M_sexChr_X_testP_df <- data.frame(M_sexChr_X_testP)
names(M_sexChr_X_testP_df)[names(M_sexChr_X_testP_df) == "M_sexChr_X_testP"] <- "pvalue"
M_sexChr_Y_testP_df <- data.frame(M_sexChr_Y_testP)
names(M_sexChr_Y_testP_df)[names(M_sexChr_Y_testP_df) == "M_sexChr_Y_testP"] <- "pvalue"



corr_allgenes_pvalues <- rbind(M_testP_df, M_sexChr_X_testP_df, M_sexChr_Y_testP_df, M_autosomal_all_testP_df)

# multiple test corrections 
M_sexdiff <- as.numeric(unlist(M_testP_df))
M_sexdiff.p.adj <- p.adjust(M_sexdiff, method = "hochberg")
M_sexdiff.p.adj <- data.frame(M_sexdiff.p.adj)
names(M_sexdiff.p.adj)[names(M_sexdiff.p.adj) == "M_sexdiff.p.adj"] <- "padjvalue"

# M_sexChr <- as.numeric(unlist(M_sexChr_testP_df))
# M_sexChr.p.adj <- p.adjust(M_sexChr, method = "hochberg")
# M_sexChr.p.adj <- data.frame(M_sexChr.p.adj)
# names(M_sexChr.p.adj)[names(M_sexChr.p.adj) == "M_sexChr.p.adj"] <- "padjvalue"

M_sexChr_X <- as.numeric(unlist(M_sexChr_X_testP_df))
M_sexChr_X.p.adj <- p.adjust(M_sexChr_X, method = "hochberg")
M_sexChr_X.p.adj <- data.frame(M_sexChr_X.p.adj)
names(M_sexChr_X.p.adj)[names(M_sexChr_X.p.adj) == "M_sexChr_X.p.adj"] <- "padjvalue"

M_sexChr_Y <- as.numeric(unlist(M_sexChr_Y_testP_df))
M_sexChr_Y.p.adj <- p.adjust(M_sexChr_Y, method = "hochberg")
M_sexChr_Y.p.adj <- data.frame(M_sexChr_Y.p.adj)
names(M_sexChr_Y.p.adj)[names(M_sexChr_Y.p.adj) == "M_sexChr_Y.p.adj"] <- "padjvalue"


M_autosomal <- as.numeric(unlist(M_autosomal_all_testP_df))
M_autosomal.p.adj <- p.adjust(M_autosomal, method = "hochberg")
M_autosomal.p.adj <- data.frame(M_autosomal.p.adj)
names(M_autosomal.p.adj)[names(M_autosomal.p.adj) == "M_autosomal.p.adj"] <- "padjvalue"

corr_allgenes_padjvalues <- rbind(M_sexdiff.p.adj, M_sexChr_X.p.adj, M_sexChr_Y.p.adj, M_autosomal.p.adj)

termCorr_pvalue <- cbind(termCorr, corr_allgenes_pvalues)
termCorr_padjvalue <- cbind(termCorr, corr_allgenes_padjvalues)

#Turn your 'treatment' column into a character vector

#termCorr_pvalue <- termCorr_pvalue[order(termCorr_pvalue$Tissues),]
#Then turn it back into a factor with the levels in the correct order
#levels <- as.character(termCorr_pvalue$Tissues)
#targets <- tibble(Tissues = c("Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain Aamygdala", "Brain Anterior cingulate cortex", "Brain Ccaudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex BA9", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood", "Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain Aamygdala", "Brain Anterior cingulate cortex", "Brain Ccaudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex BA9", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood", "Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain Aamygdala", "Brain Anterior cingulate cortex", "Brain Ccaudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex BA9", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood"))
targets <- tibble(Tissues = c("Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain amygdala", "Brain anterior cingulate cortex", "Brain caudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood",
                              "Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain amygdala", "Brain anterior cingulate cortex", "Brain caudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood",
                              "Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain amygdala", "Brain anterior cingulate cortex", "Brain caudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood", 
                              "Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain amygdala", "Brain anterior cingulate cortex", "Brain caudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood"))

#df <- termCorr_pvalue[match(levels, termCorr_pvalue$Tissues),]

df <- left_join(data.frame(name=targets),termCorr_pvalue,by="Tissues")
df$Tissues <- factor(df$Tissues, levels = rev(unique(df$Tissues)))
PanelA_all_genes_df <- subset(df, df$corr == "All genes")


#------
p <- ggplot(df, aes(x = corr, y = Tissues, fill = value)) +
  scale_fill_gradient2(low = "red", mid="white", high = "steelblue", limits = c(-1, 1.00)) +
  geom_tile() + 
  ylab("")
p
r2Plot <- p +  geom_text(aes(label=value)) +  theme(axis.title.x=element_blank(),
                                                    axis.text.x=element_text(angle = 45, hjust =1),
                                                    axis.ticks.x=element_blank())

r2Plot_allGenes <- r2Plot + geom_text(aes(label= value), fontface = ifelse(df$pvalue < 0.05, 2, 1), colour = ifelse(df$pvalue < 0.05, "black", "gray40"))
r2Plot_allGenes



ggsave("PanelA_r2_allGenes.pdf", 
       r2Plot_allGenes,
       width = 5.5, height = 9.5)
dev.off()
dev.off()


#-- with multiple test corrections 

targets <- tibble(Tissues = c("Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain amygdala", "Brain anterior cingulate cortex", "Brain caudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood",
                              "Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain amygdala", "Brain anterior cingulate cortex", "Brain caudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood",
                              "Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain amygdala", "Brain anterior cingulate cortex", "Brain caudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood"))

df <- left_join(data.frame(name=targets),termCorr_padjvalue,by="Tissues")
df$Tissues <- factor(df$Tissues, levels = rev(unique(df$Tissues)))
PanelA_all_genes_df <- subset(df, df$corr == "All genes")

p <- ggplot(df, aes(x = corr, y = Tissues, fill = value)) +
  scale_fill_gradient2(low = "red", mid="white", high = "steelblue", limits = c(-1, 1.00)) +
  geom_tile() + 
  ylab("")
p
r2Plot <- p +  geom_text(aes(label=value)) +  theme(axis.title.x=element_blank(),
                                                    axis.text.x=element_text(angle = 45, hjust =1),
                                                    axis.ticks.x=element_blank())

r2Plot_allGenes <- r2Plot + geom_text(aes(label= value), fontface = ifelse(df$padjvalue < 0.05, 2, 1), colour = ifelse(df$padjvalue < 0.05, "black", "gray40"))
r2Plot_allGenes

ggsave("PanelA_r2_allGenes_multipleTest.pdf", 
       r2Plot_allGenes,
       width = 5.5, height = 9.5)
dev.off()
dev.off()

#------ EXCLUDING PLACENTA SEX DIFFERNTIALLY EXPRESSED GENES
#-----
# all genes, sex chromosome linked genes, autosomal genes
#-----
tissueMatrix <- as.matrix(df_tissues_noNA)
# dim(tissueMatrix)
sexDiffGenes.term <- subset(term_placenta, term_placenta$adj.P.Val < 0.05)
sexDiffGenes.first <- subset(firstTrimester_placenta, firstTrimester_placenta$adj.P.Val < 0.05)
sexDiffGenes <- rbind(sexDiffGenes.term, sexDiffGenes.first)
tissueMatrix <- subset(tissueMatrix, !(row.names(tissueMatrix) %in% sexDiffGenes$Geneid))
M <- cor(tissueMatrix)
M_test <- rcorr(tissueMatrix)

M_testP <- as.data.frame(M_test$P)
M_testP <- as.numeric(M_testP$`Term placenta`)
M_df <- as.data.frame(M)
M_df$`Term placenta` <- as.numeric(M_df$`Term placenta`)
Mdf_round <- M_df %>% mutate_if(is.numeric, round, digits=3)
Mdf_round$Tissues <- row.names(M_df)
Mdf_round$Tissues <- factor(Mdf_round$Tissues, levels = Mdf_round$Tissues[order(Mdf_round$`Term placenta`)])

dim(tissueMatrix)
tissueMatrix <- as.matrix(df_tissues_noNA)

# chr X
sexDiffGenes_X <- subset(sexDiffGenes, chr == "chrX")
sexChrGenes_X <- subset(term_placenta, chr == "chrX")
X_linked_notSD <- subset(sexChrGenes_X, !(sexChrGenes_X$Geneid %in% sexDiffGenes_X$Geneid))
sexChrGenes_X_tissueMatrix <- subset(tissueMatrix, (row.names(tissueMatrix) %in% X_linked_notSD$Geneid))
dim(sexChrGenes_X_tissueMatrix)

M_sexChr_X <- cor(sexChrGenes_X_tissueMatrix)
M_sexChr_X_test <- rcorr(sexChrGenes_X_tissueMatrix)
#M_sexChr_test$P
M_sexChr_X_testP <- as.data.frame(M_sexChr_X_test$P)
M_sexChr_X_testP <- as.numeric(M_sexChr_X_testP$`Term placenta`)
M_sexChr_X_df <- as.data.frame(M_sexChr_X)
M_sexChr_X_df$`Term placenta` <- as.numeric(M_sexChr_X_df$`Term placenta`)
M_sexChr_X_df_round <- M_sexChr_X_df %>% mutate_if(is.numeric, round, digits=3)
M_sexChr_X_df_round$Tissues <- row.names(M_sexChr_X_df)
M_sexChr_X_df_round$Tissues <- factor(M_sexChr_X_df_round$Tissues, levels = M_sexChr_X_df_round$Tissues[order(M_sexChr_X_df_round$`Term placenta`)])

# chr Y
sexDiffGenes_Y <- subset(sexDiffGenes, chr == "chrY")
sexChrGenes_Y <- subset(term_placenta, chr == "chrY")
Y_linked_notSD <- subset(sexChrGenes_Y, !(sexChrGenes_Y$Geneid %in% sexDiffGenes_Y$Geneid))
sexChrGenes_Y_tissueMatrix <- subset(tissueMatrix, (row.names(tissueMatrix) %in% Y_linked_notSD$Geneid))
dim(sexChrGenes_Y_tissueMatrix)
M_sexChr_Y <- cor(sexChrGenes_Y_tissueMatrix)
M_sexChr_Y_test <- rcorr(sexChrGenes_Y_tissueMatrix)
#M_sexChr_test$P
M_sexChr_Y_testP <- as.data.frame(M_sexChr_Y_test$P)
M_sexChr_Y_testP <- as.numeric(M_sexChr_Y_testP$`Term placenta`)
M_sexChr_Y_df <- as.data.frame(M_sexChr_Y)
M_sexChr_Y_df$`Term placenta` <- as.numeric(M_sexChr_Y_df$`Term placenta`)
M_sexChr_Y_df_round <- M_sexChr_Y_df %>% mutate_if(is.numeric, round, digits=3)
M_sexChr_Y_df_round$Tissues <- row.names(M_sexChr_Y_df)
M_sexChr_Y_df_round$Tissues <- factor(M_sexChr_Y_df_round$Tissues, levels = M_sexChr_Y_df_round$Tissues[order(M_sexChr_Y_df_round$`Term placenta`)])

#-------------
#tissueMatrix <- as.matrix(df_tissues_noNA)
autosomal_all <- subset(term_placenta, chr != "chrX" & chr != "chrY")
tissueMatrix <- subset(tissueMatrix, !(row.names(tissueMatrix) %in% sexDiffGenes$Geneid))
autosomal_all_tissueMatrix <- subset(tissueMatrix, (row.names(tissueMatrix) %in% autosomal_all$Geneid))
dim(autosomal_all_tissueMatrix)
M_autosomal_all <- cor(autosomal_all_tissueMatrix)
M_autosomal_all_test <- rcorr(M_autosomal_all)
#M_autosomal_all_test$P
M_autosomal_all_testP <- as.data.frame(M_autosomal_all_test$P)
M_autosomal_all_testP_rownames <- rownames(M_autosomal_all_testP)
M_autosomal_all_testP <- as.numeric(M_autosomal_all_testP$`Term placenta`)
df_M_autorownames <- as.data.frame(M_autosomal_all_testP_rownames)
df_M_autosomal_all_testP <- as.data.frame(M_autosomal_all_testP)
rownames(df_M_autosomal_all_testP) <- df_M_autorownames$M_autosomal_all_testP_rownames
M_autosomal_all_df <- as.data.frame(M_autosomal_all)
M_autosomal_all_df$`Term placenta` <- as.numeric(M_autosomal_all_df$`Term placenta`)
M_autosomal_all_df_round <- M_autosomal_all_df %>% mutate_if(is.numeric, round, digits=3)
M_autosomal_all_df_round$Tissues <- row.names(M_autosomal_all_df)
M_autosomal_all_df_round$Tissues <- factor(M_autosomal_all_df_round$Tissues, levels = M_autosomal_all_df_round$Tissues[order(M_autosomal_all_df_round$`Term placenta`)])

#--- all genes, sex-chromosome linked and autosomal 
Mdf_round$corr <- "All genes"
#M_sexChr_df_round$corr <- "Sex-linked"
M_autosomal_all_df_round$corr <- "Autosomal"
M_sexChr_X_df_round$corr <- "X-linked"
M_sexChr_Y_df_round$corr <- "Y-linked"


groupCorr <- rbind(Mdf_round, M_sexChr_X_df_round, M_sexChr_Y_df_round, M_autosomal_all_df_round)
library(reshape2)
meltCorr <- melt(groupCorr)
termCorr <- subset(meltCorr, meltCorr$variable == "Term placenta")
termCorr$pvalue <- NULL

# Pvalues 
M_testP_df <- data.frame(M_testP)
names(M_testP_df)[names(M_testP_df) == "M_testP"] <- "pvalue"
# M_sexChr_testP_df <- data.frame(M_sexChr_testP)
# names(M_sexChr_testP_df)[names(M_sexChr_testP_df) == "M_sexChr_testP"] <- "pvalue"
M_autosomal_all_testP_df <- data.frame(M_autosomal_all_testP)
names(M_autosomal_all_testP_df)[names(M_autosomal_all_testP_df) == "M_autosomal_all_testP"] <- "pvalue"
M_sexChr_X_testP_df <- data.frame(M_sexChr_X_testP)
names(M_sexChr_X_testP_df)[names(M_sexChr_X_testP_df) == "M_sexChr_X_testP"] <- "pvalue"
M_sexChr_Y_testP_df <- data.frame(M_sexChr_Y_testP)
names(M_sexChr_Y_testP_df)[names(M_sexChr_Y_testP_df) == "M_sexChr_Y_testP"] <- "pvalue"



corr_allgenes_pvalues <- rbind(M_testP_df, M_sexChr_X_testP_df, M_sexChr_Y_testP_df, M_autosomal_all_testP_df)

# multiple test corrections 
M_sexdiff <- as.numeric(unlist(M_testP_df))
M_sexdiff.p.adj <- p.adjust(M_sexdiff, method = "hochberg")
M_sexdiff.p.adj <- data.frame(M_sexdiff.p.adj)
names(M_sexdiff.p.adj)[names(M_sexdiff.p.adj) == "M_sexdiff.p.adj"] <- "padjvalue"

# M_sexChr <- as.numeric(unlist(M_sexChr_testP_df))
# M_sexChr.p.adj <- p.adjust(M_sexChr, method = "hochberg")
# M_sexChr.p.adj <- data.frame(M_sexChr.p.adj)
# names(M_sexChr.p.adj)[names(M_sexChr.p.adj) == "M_sexChr.p.adj"] <- "padjvalue"

M_sexChr_X <- as.numeric(unlist(M_sexChr_X_testP_df))
M_sexChr_X.p.adj <- p.adjust(M_sexChr_X, method = "hochberg")
M_sexChr_X.p.adj <- data.frame(M_sexChr_X.p.adj)
names(M_sexChr_X.p.adj)[names(M_sexChr_X.p.adj) == "M_sexChr_X.p.adj"] <- "padjvalue"

M_sexChr_Y <- as.numeric(unlist(M_sexChr_Y_testP_df))
M_sexChr_Y.p.adj <- p.adjust(M_sexChr_Y, method = "hochberg")
M_sexChr_Y.p.adj <- data.frame(M_sexChr_Y.p.adj)
names(M_sexChr_Y.p.adj)[names(M_sexChr_Y.p.adj) == "M_sexChr_Y.p.adj"] <- "padjvalue"


M_autosomal <- as.numeric(unlist(M_autosomal_all_testP_df))
M_autosomal.p.adj <- p.adjust(M_autosomal, method = "hochberg")
M_autosomal.p.adj <- data.frame(M_autosomal.p.adj)
names(M_autosomal.p.adj)[names(M_autosomal.p.adj) == "M_autosomal.p.adj"] <- "padjvalue"

corr_allgenes_padjvalues <- rbind(M_sexdiff.p.adj, M_sexChr_X.p.adj, M_sexChr_Y.p.adj, M_autosomal.p.adj)

termCorr_pvalue <- cbind(termCorr, corr_allgenes_pvalues)
termCorr_padjvalue <- cbind(termCorr, corr_allgenes_padjvalues)

#Turn your 'treatment' column into a character vector

#termCorr_pvalue <- termCorr_pvalue[order(termCorr_pvalue$Tissues),]
#Then turn it back into a factor with the levels in the correct order
#levels <- as.character(termCorr_pvalue$Tissues)
#targets <- tibble(Tissues = c("Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain Aamygdala", "Brain Anterior cingulate cortex", "Brain Ccaudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex BA9", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood", "Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain Aamygdala", "Brain Anterior cingulate cortex", "Brain Ccaudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex BA9", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood", "Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain Aamygdala", "Brain Anterior cingulate cortex", "Brain Ccaudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex BA9", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood"))
targets <- tibble(Tissues = c("Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain amygdala", "Brain anterior cingulate cortex", "Brain caudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood",
                              "Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain amygdala", "Brain anterior cingulate cortex", "Brain caudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood",
                              "Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain amygdala", "Brain anterior cingulate cortex", "Brain caudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood", 
                              "Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain amygdala", "Brain anterior cingulate cortex", "Brain caudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood"))

#df <- termCorr_pvalue[match(levels, termCorr_pvalue$Tissues),]

df <- left_join(data.frame(name=targets),termCorr_pvalue,by="Tissues")
df$Tissues <- factor(df$Tissues, levels = rev(unique(df$Tissues)))



#------
p <- ggplot(df, aes(x = corr, y = Tissues, fill = value)) +
  scale_fill_gradient2(low = "red", mid="white", high = "steelblue", limits = c(-1, 1.00)) +
  geom_tile() + 
  ylab("")
p
r2Plot <- p +  geom_text(aes(label=value)) +  theme(axis.title.x=element_blank(),
                                                    axis.text.x=element_text(angle = 45, hjust =1),
                                                    axis.ticks.x=element_blank())

r2Plot_allGenes <- r2Plot + geom_text(aes(label= value), fontface = ifelse(df$pvalue < 0.05, 2, 1), colour = ifelse(df$pvalue < 0.05, "black", "gray40"))
r2Plot_allGenes



ggsave("PanelB_r2_allGenes_excludingPlacentaSDgenes.pdf", 
       r2Plot_allGenes,
       width = 5.5, height = 9.5)
dev.off()
dev.off()


#-- with multiple test corrections 

targets <- tibble(Tissues = c("Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain amygdala", "Brain anterior cingulate cortex", "Brain caudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood",
                              "Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain amygdala", "Brain anterior cingulate cortex", "Brain caudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood",
                              "Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain amygdala", "Brain anterior cingulate cortex", "Brain caudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood"))

df <- left_join(data.frame(name=targets),termCorr_padjvalue,by="Tissues")
df$Tissues <- factor(df$Tissues, levels = rev(unique(df$Tissues)))
PanelB_all_genes_df <- subset(df, df$corr == "All genes")

p <- ggplot(df, aes(x = corr, y = Tissues, fill = value)) +
  scale_fill_gradient2(low = "red", mid="white", high = "steelblue", limits = c(-1, 1.00)) +
  geom_tile() + 
  ylab("")
p
r2Plot <- p +  geom_text(aes(label=value)) +  theme(axis.title.x=element_blank(),
                                                    axis.text.x=element_text(angle = 45, hjust =1),
                                                    axis.ticks.x=element_blank())

r2Plot_allGenes <- r2Plot + geom_text(aes(label= value), fontface = ifelse(df$padjvalue < 0.05, 2, 1), colour = ifelse(df$padjvalue < 0.05, "black", "gray40"))
r2Plot_allGenes

ggsave("PanelB_r2_allGenes_multipleTest_excludingPlacentaSDgenes.pdf", 
       r2Plot_allGenes,
       width = 5.5, height = 9.5)
dev.off()
dev.off()



#-----
# placenta sex deferentially expressed genes, sex chromosome linked genes, autosomal genes
#-----
# corr of sex deferentially expressed genes 
sexDiffGenes <- subset(term_placenta, term_placenta$adj.P.Val < 0.05)
sexDiffGenesFirstTri <- subset(firstTrimester_placenta, firstTrimester_placenta$adj.P.Val < 0.05)
sexDiffGenes <- rbind(sexDiffGenes, sexDiffGenesFirstTri)
sexDiffGeneCount <- unique(sexDiffGenes$Geneid)
df_tissues_noNA <- na.omit(df_tissues)
dim(df_tissues)
subset(df_tissues, rownames(df_tissues) %in% "MANEA")
df_tissues_noNA_termSexDiff <- subset(df_tissues_noNA, row.names(df_tissues_noNA) %in% sexDiffGenes$Geneid)
# correlation
tissueMatrix_termSexDiff <- as.matrix(df_tissues_noNA_termSexDiff)
dim(tissueMatrix_termSexDiff)
# All sex differentially expressed genes 
M_SD <- cor(tissueMatrix_termSexDiff)
M_SD_df <- as.data.frame(M_SD)
M_SD_rcorr <- rcorr(tissueMatrix_termSexDiff)
M_SD_rcorr_P <- as.data.frame(M_SD_rcorr$P)
M_SD_rcorr_P$`Term placenta` <- as.numeric(M_SD_rcorr_P$`Term placenta`)
M_SD_df$`Term placenta` <- as.numeric(M_SD_df$`Term placenta`)
M_SD_df_round <- M_SD_df %>% mutate_if(is.numeric, round, digits=3)
M_SD_df_round$Tissues <- row.names(M_SD_df)
M_SD_df_round$Tissues <- factor(M_SD_df_round$Tissues, levels = M_SD_df_round$Tissues[order(M_SD_df_round$`Term placenta`)])

# sex linked genes 
sexChrGenes_SD_X <- subset(sexDiffGenes, chr == "chrX")
sexChrGenes_SD_X_tissueMatrix <- subset(tissueMatrix_termSexDiff, (row.names(tissueMatrix_termSexDiff) %in% sexChrGenes_SD_X$Geneid))
dim(sexChrGenes_SD_X_tissueMatrix)
M_sexChr_SD_X <- cor(sexChrGenes_SD_X_tissueMatrix)
M_sexChr_df_X <- as.data.frame(M_sexChr_SD_X)
M_sexChr_df_X$`Term placenta` <- as.numeric(M_sexChr_df_X$`Term placenta`)
M_sexChr_X_rcorr <- rcorr(sexChrGenes_SD_X_tissueMatrix)
M_sexChr_X_rcorr_P <- as.data.frame(M_SD_rcorr$P)
M_sexChr_X_rcorr_P$`Term placenta` <- as.numeric(M_sexChr_X_rcorr_P$`Term placenta`)
M_sexChr_X_df_round <- M_sexChr_df_X %>% mutate_if(is.numeric, round, digits=3)
M_sexChr_X_df_round$Tissues <- row.names(M_sexChr_df_X)
M_sexChr_X_df_round$Tissues <- factor(M_sexChr_X_df_round$Tissues, levels = M_sexChr_X_df_round$Tissues[order(M_sexChr_X_df_round$`Term placenta`)])
dim(sexChrGenes_SD_X_tissueMatrix)

# chr Y
sexChrGenes_SD_Y <- subset(sexDiffGenes, chr == "chrY")
sexChrGenes_SD_Y_tissueMatrix <- subset(tissueMatrix_termSexDiff, (row.names(tissueMatrix_termSexDiff) %in% sexChrGenes_SD_Y$Geneid))
dim(sexChrGenes_SD_Y_tissueMatrix)
M_sexChr_SD_Y <- cor(sexChrGenes_SD_Y_tissueMatrix)
M_sexChr_df_Y <- as.data.frame(M_sexChr_SD_Y)
M_sexChr_df_Y$`Term placenta` <- as.numeric(M_sexChr_df_Y$`Term placenta`)
M_sexChr_Y_rcorr <- rcorr(sexChrGenes_SD_Y_tissueMatrix)
M_sexChr_Y_rcorr_P <- as.data.frame(M_SD_rcorr$P)
M_sexChr_Y_rcorr_P$`Term placenta` <- as.numeric(M_sexChr_Y_rcorr_P$`Term placenta`)
M_sexChr_Y_df_round <- M_sexChr_df_Y %>% mutate_if(is.numeric, round, digits=3)
M_sexChr_Y_df_round$Tissues <- row.names(M_sexChr_df_Y)
M_sexChr_Y_df_round$Tissues <- factor(M_sexChr_Y_df_round$Tissues, levels = M_sexChr_Y_df_round$Tissues[order(M_sexChr_Y_df_round$`Term placenta`)])
dim(sexChrGenes_SD_Y_tissueMatrix)
# autosomal genes 
autosomal_all_SD <- subset(sexDiffGenes, chr != "chrX" & chr != "chrY")
autosomal_all_SD_tissueMatrix <- subset(tissueMatrix_termSexDiff, (row.names(tissueMatrix_termSexDiff) %in% autosomal_all_SD$Geneid))
dim(autosomal_all_SD_tissueMatrix)
M_autosomal_all_SD <- cor(autosomal_all_SD_tissueMatrix)
M_autosomal_all_SD_df <- as.data.frame(M_autosomal_all_SD)
M_autosomal_all_SD_df$`Term placenta` <- as.numeric(M_autosomal_all_SD_df$`Term placenta`)
M_autosomal_rcorr <- rcorr(autosomal_all_SD_tissueMatrix)
M_autosomal_rcorr_P <- as.data.frame(M_SD_rcorr$P)
M_autosomal_rcorr_P$`Term placenta` <- as.numeric(M_autosomal_rcorr_P$`Term placenta`)
dim(autosomal_all_SD_tissueMatrix)

M_autosomal_all_SD_df_round <- M_autosomal_all_SD_df %>% mutate_if(is.numeric, round, digits=3)
M_autosomal_all_SD_df_round$Tissues <- row.names(M_autosomal_all_SD_df)
M_autosomal_all_SD_df_round$Tissues <- factor(M_autosomal_all_SD_df_round$Tissues, levels = M_autosomal_all_SD_df_round$Tissues[order(M_autosomal_all_SD_df_round$`Term placenta`)])
#--- all genes, sex-chromosome linked and autosomal 
M_SD_df_round$corr <- "All genes"
M_sexChr_X_df_round$corr <- "X-linked"
M_sexChr_Y_df_round$corr <- "Y-linked"

M_autosomal_all_SD_df_round$corr <- "Autosomal"

groupCorr <- rbind(M_SD_df_round, M_sexChr_X_df_round, M_sexChr_Y_df_round, M_autosomal_all_SD_df_round)
library(reshape2)
meltCorr <- melt(groupCorr)
termCorr <- subset(meltCorr, meltCorr$variable == "Term placenta")

M_SD_rcorr_P <- data.frame(M_SD_rcorr_P$`Term placenta`)
names(M_SD_rcorr_P)[names(M_SD_rcorr_P) == "M_SD_rcorr_P..Term.placenta."] <- "pvalue"
M_sexChr_X_rcorr_P <- data.frame(M_sexChr_X_rcorr_P$`Term placenta`)
names(M_sexChr_X_rcorr_P)[names(M_sexChr_X_rcorr_P) == "M_sexChr_X_rcorr_P..Term.placenta."] <- "pvalue"
M_sexChr_Y_rcorr_P <- data.frame(M_sexChr_Y_rcorr_P$`Term placenta`)
names(M_sexChr_Y_rcorr_P)[names(M_sexChr_Y_rcorr_P) == "M_sexChr_Y_rcorr_P..Term.placenta."] <- "pvalue"
M_autosomal_rcorr_P <- data.frame(M_autosomal_rcorr_P$`Term placenta`)
names(M_autosomal_rcorr_P)[names(M_autosomal_rcorr_P) == "M_autosomal_rcorr_P..Term.placenta."] <- "pvalue"

corr_SDgenes_pvalues <- rbind(M_SD_rcorr_P, M_sexChr_X_rcorr_P, M_sexChr_Y_rcorr_P, M_autosomal_rcorr_P)
termCorr_pvalue <- cbind(termCorr, corr_SDgenes_pvalues)

#---# multiple test corrections 
# multiple test corrections 
M_sexdiff <- as.numeric(unlist(M_SD_rcorr_P))
M_sexdiff.p.adj <- p.adjust(M_sexdiff, method = "hochberg")
M_sexdiff.p.adj <- data.frame(M_sexdiff.p.adj)
names(M_sexdiff.p.adj)[names(M_sexdiff.p.adj) == "M_sexdiff.p.adj"] <- "padjvalue"

M_sexChr_X <- as.numeric(unlist(M_sexChr_X_rcorr_P))
M_sexChr_X.p.adj <- p.adjust(M_sexChr_X, method = "hochberg")
M_sexChr_X.p.adj <- data.frame(M_sexChr_X.p.adj)
names(M_sexChr_X.p.adj)[names(M_sexChr_X.p.adj) == "M_sexChr_X.p.adj"] <- "padjvalue"

M_sexChr_Y <- as.numeric(unlist(M_sexChr_Y_rcorr_P))
M_sexChr_Y.p.adj <- p.adjust(M_sexChr_Y, method = "hochberg")
M_sexChr_Y.p.adj <- data.frame(M_sexChr_Y.p.adj)
names(M_sexChr_Y.p.adj)[names(M_sexChr_Y.p.adj) == "M_sexChr_Y.p.adj"] <- "padjvalue"

M_autosomal <- as.numeric(unlist(M_autosomal_rcorr_P))
M_autosomal.p.adj <- p.adjust(M_autosomal, method = "hochberg")
M_autosomal.p.adj <- data.frame(M_autosomal.p.adj)
names(M_autosomal.p.adj)[names(M_autosomal.p.adj) == "M_autosomal.p.adj"] <- "padjvalue"

corr_allgenes_padjvalues <- rbind(M_sexdiff.p.adj, M_sexChr_X.p.adj, M_sexChr_Y.p.adj, M_autosomal.p.adj)
termCorr_padjvalue <- cbind(termCorr, corr_allgenes_padjvalues)


#Turn your 'treatment' column into a character vector

#termCorr_pvalue <- termCorr_pvalue[order(termCorr_pvalue$Tissues),]
#Then turn it back into a factor with the levels in the correct order
#levels <- as.character(termCorr_pvalue$Tissues)
#targets <- tibble(Tissues = c("Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain amygdala", "Brain anterior cingulate cortex", "Brain caudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood", "Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain amygdala", "Brain Anterior cingulate cortex", "Brain Ccaudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex BA9", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood", "Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain Aamygdala", "Brain Anterior cingulate cortex", "Brain Ccaudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex BA9", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood"))
targets <- tibble(Tissues = c("Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain amygdala", "Brain anterior cingulate cortex", "Brain caudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood",
"Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain amygdala", "Brain anterior cingulate cortex", "Brain caudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood",
"Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain amygdala", "Brain anterior cingulate cortex", "Brain caudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood",
"Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain amygdala", "Brain anterior cingulate cortex", "Brain caudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood"))

#df <- termCorr_pvalue[match(levels, termCorr_pvalue$Tissues),]

df <- left_join(data.frame(name=targets),termCorr_pvalue,by="Tissues")
df$Tissues <- factor(df$Tissues, levels = rev(unique(df$Tissues)))

#------
p <- ggplot(df, aes(x = corr, y = Tissues, fill = value)) +
  scale_fill_gradient2(low = "red", mid="white", high = "steelblue", limits = c(-1, 1.00)) +
  geom_tile() + 
  ylab("")
p
r2Plot <- p +  geom_text(aes(label=value)) +  theme(axis.title.x=element_blank(),
                                                    axis.text.x=element_text(angle = 45, hjust =1),
                                                    axis.ticks.x=element_blank())

r2Plot_allGenes <- r2Plot + geom_text(aes(label= value), fontface = ifelse(df$pvalue < 0.05, 2, 1), colour = ifelse(df$pvalue < 0.05, "black", "gray40"))
ggsave("PanelC_r2_SDGenes.pdf", 
       r2Plot_allGenes,
       width = 5.5, height = 9.5)
dev.off()
dev.off()

targets <- tibble(Tissues = c("Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain amygdala", "Brain anterior cingulate cortex", "Brain caudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood",
                              "Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain amygdala", "Brain anterior cingulate cortex", "Brain caudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood",
                              "Term placenta", "Late first trimester placenta", "Adipose subcutaneous", "Adipose visceral omentum", "Adrenal gland", "Artery aorta", "Artery coronary", "Artery tibial", "Bladder", "Brain amygdala", "Brain anterior cingulate cortex", "Brain caudate", "Brain cerebellar hemisphere", "Brain cerebellum", "Brain cortex", "Brain frontal cortex", "Brain hippocampus", "Brain hypothalamus", "Brain nucleus accumbens", "Brain putamen", "Brain spinal cord", "Brain substantia nigra", "Breast mammary", "Colon", "Colon transverse", "Esophagus", "Esophagus mucosa", "Esophagus muscularis", "Heart atrial appendage", "Heart left ventricle", "Kidney", "Liver", "Lung", "Minor salivary gland", "Nerve tibial", "Pancreas", "Pituitary", "Skin not sun exposed suprapubic", "Skin sun exposed lower leg", "Small intestine", "Spleen", "Stomach", "Thyroid", "Whole blood"))

df <- left_join(data.frame(name=targets),termCorr_padjvalue,by="Tissues")
df$Tissues <- factor(df$Tissues, levels = rev(unique(df$Tissues)))

p <- ggplot(df, aes(x = corr, y = Tissues, fill = value)) +
  scale_fill_gradient2(low = "red", mid="white", high = "steelblue", limits = c(-1, 1.00)) +
  geom_tile() + 
  ylab("")
p
r2Plot <- p +  geom_text(aes(label=value)) +  theme(axis.title.x=element_blank(),
                                                    axis.text.x=element_text(angle = 45, hjust =1),
                                                    axis.ticks.x=element_blank())

r2Plot_allGenes <- r2Plot + geom_text(aes(label= value), fontface = ifelse(df$padjvalue < 0.05, 2, 1), colour = ifelse(df$padjvalue < 0.05, "black", "gray40"))
r2Plot_allGenes

ggsave("PanelC_r2_SDGenes_multipleTest.pdf", 
       r2Plot_allGenes,
       width = 5.5, height = 9.5)
dev.off()
dev.off()

#-------------- main figure plot

PanelB_all_genes_df$corr <- gsub('All genes', 'Exclude SD genes', PanelB_all_genes_df$corr)
df$corr <- gsub('All genes', 'All SD genes', df$corr)
main_fig <- rbind(PanelA_all_genes_df, PanelB_all_genes_df, df)

#df$Tissues <- factor(df$Tissues, levels = rev(unique(df$Tissues)))

unique(main_fig$corr)
main_fig$corr <- factor(main_fig$corr,levels = c("All genes", "Exclude SD genes", "All SD genes", "Autosomal", "X-linked", "Y-linked"))
#------
p <- ggplot(main_fig, aes(x = corr, y = Tissues, fill = value)) +
  scale_fill_gradient2(low = "red", mid="white", high = "steelblue", limits = c(-1, 1.00)) +
  geom_tile() + 
  ylab("")
p
r2Plot <- p +  geom_text(aes(label=sprintf("%0.3f", round(value, digits = 3)))) +  theme(axis.title.x=element_blank(),
                                                    axis.text.x=element_text(angle = 45, hjust =1),
                                                    axis.ticks.x=element_blank())
r2Plot
r2Plot_allGenes <- r2Plot + geom_text(aes(label=sprintf("%0.3f", round(value, digits = 3))), fontface = ifelse(main_fig$padjvalue < 0.05, 2, 1), colour = ifelse(main_fig$padjvalue < 0.05, "black", "gray40"))
r2Plot_allGenes <- r2Plot_allGenes + scale_x_discrete(labels=c("All genes" = "All", "Exclude SD genes" = "All exclude SD",  "All SD genes" = "All SD", "Autosomal" = "Autosomal SD", "X-linked" = "X-linked SD", "Y-linked" = "Y-linked SD"))
#r2Plot_allGenes + geom_rect(aes(xmin = "All genes", xmax = "Y-linked", ymin = "Brain amygdala", ymax = "Brain substantia nigra"), fill = "transparent", color = "red", size = 1.5)

r2Plot_allGenes
# convert x and y variables to factors
corr <- as.factor(main_fig$corr)
Tissues <- as.factor(main_fig$Tissues)
# numeric version of the levels to be bound by a box
xmin <- unique(as.numeric(corr[corr == "All genes"]))
xmax <- unique(as.numeric(corr[corr == "Y-linked"]))

ymin <- unique(as.numeric(Tissues[Tissues == "Bladder"]))
ymax <- unique(as.numeric(Tissues[Tissues == "Breast mammary"]))

# set offset
offset <- 0.5

new <- r2Plot_allGenes + geom_rect(aes(xmin = xmin - offset,
              xmax = xmax + offset,
              ymin = ymin - offset,
              ymax = ymax + offset),
          fill = "transparent", color = "black", size = .75)

new
ggsave("mainFig_r2_genes_emptySpace.pdf", 
       new,
       width = 6.5, height = 9.5)
dev.off()
dev.off()
