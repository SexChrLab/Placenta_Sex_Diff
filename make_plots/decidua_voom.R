#   title: "decidua MDS plots"
# author: "Kimberly Olney"
# date: "02/14/2020"
# output:
#   pdf_document: default
# html_document:
#   df_print: paged

library(limma)
library(edgeR)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(mixOmics)
library(devtools)
library(dplyr)
library(stringr)
library(ggrepel)
library(Glimma)
library(gplots)
library(plyr)
library(variancePartition)
library(doParallel)

# set working directory 
setwd("~/Dropbox (ASU)/Placenta/BATCH2_PLACENTA_DECIDUA_ANALYSIS/HISAT_FeatureCounts/batch1_and_batch2/")

# Read in count, genes, and phenotype data
# For gene level and transcript level analysis
genes <- read.csv("counts_pheno/genesID.csv", header=TRUE, sep = ",") # gene ids
transcripts <- read.csv("counts_pheno/transcriptsID.csv", header=TRUE, sep = ",") # transcript ids


# decidua gene analysis information (donated by 'dg')
counts_dg <- read.delim("counts_pheno/decidua_geneCounts.tsv", header=TRUE, sep="\t")
colnames(counts_dg) <- str_replace_all(colnames(counts_dg), pattern="\\.","-") # replace . with - in sample names

counts_dt <- read.delim("counts_pheno/decidua_transcriptCounts.tsv", header=TRUE, sep="\t")
colnames(counts_dt) <- str_replace_all(colnames(counts_dt), pattern="\\.","-") # replace . with - in sample names

pheno_d <- read.csv("counts_pheno/200508_decidua_pheno.csv", header=TRUE, sep=",")

# decidua samples to remove due to failed QC and/or outlier in MDS plot
decidua_removals <- c("YPOPS0007M-DEC", "OBG0021-DEC", "OBG0019-DEC", "OBG0107-DEC")

no_removals <- c()

samplesToRemove <- c(decidua_removals) # update depending on comparison being made
SAMPLE_LENGTH <- as.numeric(length(samplesToRemove)) # to call later 
half_sample_length <- SAMPLE_LENGTH/2 # half the sample length

# update counts data frame to exlude select samples 
removals_dg <- (names(counts_dg) %in% samplesToRemove[1:SAMPLE_LENGTH]) # for matching names create a value of TRUE or FALSES 
counts_dg_ExRemovals <-counts_dg[!removals_dg] # create a new counts file that excludes (Ex) the removals IDs 

removals_dt <- (names(counts_dt) %in% samplesToRemove[1:SAMPLE_LENGTH]) # for matching names create a value of TRUE or FALSES 
counts_dt_ExRemovals <-counts_dt[!removals_dt] # create a new counts file that excludes (Ex) the removals IDs 

geneCounts <- cbind(genes, counts_dg_ExRemovals)
transcriptsCounts <- cbind(transcripts, counts_dt_ExRemovals)

# update pheno data frame to exlude select samples 
pheno_ExRemovals <- pheno_d[! pheno_d$sample %in% samplesToRemove[1:SAMPLE_LENGTH],] # update 1:16 depending on size of samples to remove

# create groups for the samples  
sex <- factor(pheno_ExRemovals$sex, levels=c("female", "male"))
sample <- factor(pheno_ExRemovals$sample)
Offspring_sex <- factor(pheno_ExRemovals$Offspring_sex, levels=c("female", "male"))
tissue <- factor(pheno_ExRemovals$tissue, levels=c("decidua", "placenta"))
batch <- factor(pheno_ExRemovals$batch, levels=c("1","2"))
Maternal_age <- factor(pheno_ExRemovals$Maternal_age)
Maternal_age <- as.numeric(levels(Maternal_age))[Maternal_age]
Gestational_age <- factor(pheno_ExRemovals$Gestational_age)
Gestational_age <- as.numeric(levels(Gestational_age))[Gestational_age]
Gravidity <- factor(pheno_ExRemovals$Gravidity)
Gravidity <- as.numeric(levels(Gravidity))[Gravidity]
Parity <- factor(pheno_ExRemovals$Parity)
Parity <- as.numeric(levels(Parity))[Parity]
BirthWeight <- factor(pheno_ExRemovals$BirthWeight)
BirthWeight <- as.numeric(levels(BirthWeight))[BirthWeight]
MethodConception <- factor(pheno_ExRemovals$MethodConception, levels = c("Spontaneous", "Intrauterine_insemination", "IVF"))
Reported_race <- factor(pheno_ExRemovals$Reported_race, levels=c("Asian", "Black", "White", "Hispanic", "Other"))
site <- factor(pheno_ExRemovals$site, levels=c("RNA_1", "RNA_2"))
Lane <- factor(pheno_ExRemovals$Lane, levels=c("L002", "L003", "L004"))
Ancestry_prediction <- factor(pheno_ExRemovals$Ancestry_prediction)
PC1 <- factor(pheno_ExRemovals$PC1)
PC1 <- as.numeric(levels(PC1))[PC1]
PC2 <- factor(pheno_ExRemovals$PC2)
PC2 <- as.numeric(levels(PC2))[PC2]
Pre.Pregnancy_BMI <- factor(pheno_ExRemovals$Pre.Pregnancy_BMI)
Pre.Pregnancy_BMI <- as.numeric(levels(Pre.Pregnancy_BMI))[Pre.Pregnancy_BMI]
# create DGE list object 
DGE_gene <- DGEList(counts=counts_dg_ExRemovals, genes = genes) # create DGEList object 
samplenames <- (pheno_ExRemovals$sample) # sample names are not unique because a sample may belong to multiple groups
colnames(DGE_gene) <- samplenames

DGE_tran <- DGEList(counts=counts_dt_ExRemovals, genes = transcripts) # create DGEList object 
samplenames <- (pheno_ExRemovals$sample) # sample names are not unique because a sample may belong to multiple groups
colnames(DGE_tran) <- samplenames

# add groups to samples in dge list 
DGE_gene$samples$sex <- sex
DGE_gene$samples$sample <- sample
DGE_gene$samples$Offspring_sex <- Offspring_sex
DGE_gene$samples$Reported_race <- Reported_race
DGE_gene$samples$site <- site
DGE_gene$samples$Lane <- Lane
DGE_gene$genes <- genes
DGE_gene$samples$MethodConception <- MethodConception
DGE_gene$samples$PC1 <- PC1
DGE_gene$samples$PC2 <- PC2
DGE_gene$samples$Gravidity <- Gravidity
DGE_gene$samples$Maternal_age <- Maternal_age
DGE_gene$samples$Gestational_age <- Gestational_age
DGE_gene$samples$Parity <- Parity
DGE_gene$samples$BirthWeight <- BirthWeight 
DGE_gene$samples$Pre.Pregnancy_BMI <- Pre.Pregnancy_BMI 
# add groups to samples in dge list 
DGE_tran$samples$sex <- sex
DGE_tran$samples$sample <- sample
DGE_tran$samples$Offspring_sex <- Offspring_sex
DGE_tran$samples$Reported_race <- Reported_race
DGE_tran$samples$site <- site
DGE_tran$samples$Lane <- Lane
DGE_tran$genes <- transcripts
DGE_tran$samples$PC1 <- PC1
DGE_tran$samples$PC2 <- PC2
DGE_tran$samples$Gravidity <- Gravidity
DGE_tran$samples$Maternal_age <- Maternal_age
DGE_tran$samples$Gestational_age <- Gestational_age
DGE_tran$samples$MethodConception <- MethodConception
DGE_tran$samples$Parity <- Parity
DGE_tran$samples$BirthWeight <- BirthWeight 
DGE_tran$samples$Pre.Pregnancy_BMI <- Pre.Pregnancy_BMI 

# data pre-processing
#---- 
DGE_gene <- calcNormFactors(DGE_gene, method="TMM")

fpkm <- rpkm(DGE_gene, gene.length=DGE_gene$genes$length)
dim(fpkm)
decidua_FPKM <- cbind(genes, fpkm)
write.table(decidua_FPKM, "genelists/deciduas/decidua_geneLevel_FPKM.txt", sep = "\t", quote = FALSE, row.names = FALSE)

female_mean_fpkm <- apply(as.data.frame(fpkm)
                          [(DGE_gene$samples$Offspring_sex=="female")],
                          1, mean, na.rm=TRUE)
male_mean_fpkm <- apply(as.data.frame(fpkm)
                        [(DGE_gene$samples$Offspring_sex=="male")],
                        1, mean, na.rm=TRUE)
keep <- (female_mean_fpkm > 1.0 | male_mean_fpkm > 1.0)
table(keep)
#------ 
DGE_gene <- DGE_gene[keep,,keep.lib.sizes=FALSE]
cpm_gene <- cpm(DGE_gene) #log2-transformed counts per gene per million mapped reads (cpm).
DGE_genes <- DGE_gene$genes
may <- cbind(DGE_genes, cpm_gene)
may$length <- NULL
may2 <- melt(may, id=c("Geneid", "chr"))
pheno_p <- pheno_ExRemovals
df_merged <- merge(may2, pheno_p, by.x="variable", by.y="sample", all.x=TRUE)
write.table(df_merged, "genelists/deciduas/cpm_decidua_df_genes.txt", sep = "\t")

lcpm_gene <- cpm(cpm_gene, log=TRUE) 

keep<- rowSums(cpm_gene>0)>=0
table(keep)
dge_g_CPM <- DGE_gene[keep,, keep.lib.sizes=FALSE]
dim(dge_g_CPM)
logdge_g_CPM1 <- cpm(dge_g_CPM, log=TRUE) # log transformting the data 
dim(logdge_g_CPM1)

#---- 
DGE_tran <- calcNormFactors(DGE_tran, method="TMM")

fpkm <- rpkm(DGE_tran, gene.length=DGE_tran$genes$length)
dim(fpkm)
decidua_FPKM <- cbind(transcripts, fpkm)
write.table(decidua_FPKM, "genelists/deciduas/decidua_transLevel_FPKM.txt", sep = "\t", quote = FALSE, row.names = FALSE)

female_mean_fpkm <- apply(as.data.frame(fpkm)
                          [(DGE_tran$samples$Offspring_sex=="female")],
                          1, mean, na.rm=TRUE)
male_mean_fpkm <- apply(as.data.frame(fpkm)
                        [(DGE_tran$samples$Offspring_sex=="male")],
                        1, mean, na.rm=TRUE)
keep <- (female_mean_fpkm > 1.0 | male_mean_fpkm > 1.0)
#------ 
DGE_tran <- DGE_tran[keep,,keep.lib.sizes=FALSE]
dim(DGE_tran$genes)
cpm_tran <- cpm(DGE_tran) #log2-transformed counts per gene per million mapped reads (cpm).
lcpm_tran <- cpm(cpm_tran, log=TRUE) 

keep<- rowSums(cpm_tran>0)>=0
table(keep)
dge_t_CPM <- DGE_tran[keep,, keep.lib.sizes=FALSE]
dim(dge_t_CPM)
logdge_t_CPM1 <- cpm(dge_t_CPM, log=TRUE) # log transformting the data 
dim(logdge_t_CPM1)
#-------------------------------------------
col.offspringSex <- Offspring_sex
levels(col.offspringSex) <-  brewer.pal(nlevels(col.offspringSex), "Accent")
col.offspringSex <- as.character(col.offspringSex)

#-------------------------------------------
png(filename ="FIGURES/MDS/deciduas/deciduas_excludingFails_common_allgenes_dim1&2.png",width=5,height=5,units="in",res=1200)
plotMDS(logdge_g_CPM1, labels=Offspring_sex, gene.selection = "common", col=col.offspringSex, dim.plot = c(1,2))
title(main="All genes")
dev.off()
dev.off()

png(filename ="FIGURES/MDS/deciduas/deciduas_excludingFails_common_top100genes_dim1&2.png",width=5,height=5,units="in",res=1200)
plotMDS(logdge_g_CPM1, labels=Offspring_sex, top = 100, gene.selection = "common", col=col.offspringSex, dim.plot = c(1,3))
title(main="Top 100 genes")
dev.off()
dev.off()

png(filename ="FIGURES/MDS/deciduas/deciduas_excludingFails_common_allTranscripts_dim1&2.png",width=5,height=5,units="in",res=1200)
plotMDS(logdge_t_CPM1, labels=Offspring_sex, gene.selection = "common", col=col.offspringSex, dim.plot = c(1,2))
title(main="All transcripts")
dev.off()
dev.off()

png(filename ="FIGURES/MDS/deciduas/deciduas_excludingFails_common_top100transcripts_dim1&2.png",width=5,height=5,units="in",res=1200)
plotMDS(logdge_t_CPM1, labels=Offspring_sex, gene.selection = "common", top = 100, col=col.offspringSex, dim.plot = c(1,2))
title(main="Top 100 trasncripts")
dev.off()
dev.off()

png(filename ="FIGURES/MDS/deciduas/deciduas_excludingFails_dim1&2_PAR2x2.png",width=5,height=5,units="in",res=1200)
par(mfrow = c(2,2), mar = c(1,1,1,1) + 0.9)

plotMDS(logdge_g_CPM1, labels=Offspring_sex, gene.selection = "common", col=col.offspringSex, dim.plot = c(1,2))
title(main="All genes")

plotMDS(logdge_g_CPM1, labels=Offspring_sex, gene.selection = "common", top = 100, col=col.offspringSex, dim.plot = c(1,2))
title(main="Top 100 genes")

plotMDS(logdge_t_CPM1, labels=Offspring_sex, gene.selection = "common", col=col.offspringSex, dim.plot = c(1,2))
title(main="All transcripts")

plotMDS(logdge_t_CPM1, labels=Offspring_sex, gene.selection = "common", top = 100, col=col.offspringSex, dim.plot = c(1,2))
title(main="Top 100 transcripts")

dev.off()
dev.off()

#------------- VOOM 
#--------- gene level Voom 
group <-interaction(DGE_gene$samples$Offspring_sex) #, DGE$samples$batch, DGE$samples$Lane, DGE$samples$Reported_race) #, DGE$samples$batch, DGE$samples$Reported_race, DGE$samples$site, DGE$samples$Lane, DGE$samples)
mm <- model.matrix(~0 + group + Maternal_age + Gravidity + Lane + PC1 + PC2)
#mm <- model.matrix(~0 + group + Maternal_age + Gravidity + MethodConception + Gestational_age + Lane + PC1 + PC2)
colnames(mm) <- gsub("v\\$targets\\$Offspring_sex", "", colnames(mm))
dim(mm)

y <- voom(DGE_gene, mm, plot = T) # plot T for plotting 

# Fitting linear models in limma
fit <- lmFit(y, mm)

png(filename ="FIGURES/voomFits/deciduas/mat_Age_Gravity_Lane_PC1_PC2/decidua_voomFit.png",width=6,height=6,units="in",res=1200)
v <- voomWithQualityWeights(DGE_gene, mm, plot=TRUE)
dev.off()
dev.off()

voomCounts <- v$E
voomCountsMatrix <- data.matrix(voomCounts, rownames.force = NA)
#write.table(voomCounts, "decidua_voom_cpm_genelevel.txt", sep="\t", row.names = FALSE, quote = FALSE)
genes <- v$genes
genesVoom <- cbind(genes, voomCounts)

#-------------------------------
# variance Part
# Fit model and extract results. A linear mixed model is used
# each entry in results is a regression model fit on a single gene
form <- ~ Maternal_age + Gravidity + Parity + Pre.Pregnancy_BMI + BirthWeight + (1|MethodConception) + (1|Reported_race) + (1|Ancestry_prediction) + (1|Lane) + (1|sex)

# 2) extract variance fractions from each model fit for each gene, returns fraction of variation attributable to each variable
# Interpretation: the variance explained by each variables after correcting for all other variables

cl <- makeCluster(4)
registerDoParallel(cl)

# Fitting the linear model
varPart <- fitExtractVarPartModel(voomCounts, form, pheno_ExRemovals)

#varPart <- fitExtractVarPartModel(df, form, pheno_ExRemovals )
#write.table(varPart, "variancePart/varPart_deciduas.txt")
# to try and run this
#the brglm2 package has a method for detecting complete separation:
# install.packages("devtools")
# devtools::install_github("ikosmidis/brglm2")
# library("brglm2")
#glm(form, data = genesVoom, pheno_ExRemovals, family = binomial, method="detect_separation")
# data frame of pheno and counts....
#This should say whether complete separation occurs, and in which (combinations of) variables, e.g.
# sort variables (i.e. columns) by median fraction of variance explained
genesVP <- cbind(genes, varPart)
write.table(genesVP, "variancePart/deciduas/varPartgenes_deciduas.txt")

# violin plot of contribution of each variable to total variance
vp <- sortCols( varPart )
png(filename ="variancePart/deciduas/varPart_deciduas.png",width=6,height=6,units="in",res=1200)
plotVarPart( vp )
dev.off()
dev.off()

#---- end var part
#-------------------------------

#write.table(genes, "decidua_genes.txt", row.names = FALSE, quote = FALSE)
#v$genes
fit <- lmFit(v, mm)
png(filename ="FIGURES/voomFits/deciduas/mat_Age_Gravity_Lane_PC1_PC2/decidua_fitAmeanHist_fpkm1.png",width=6,height=6,units="in",res=1200)
hist(fit$Amean) 
dev.off()
dev.off()
plotSA(fit)
# keep <- fit$Amean > 0
# fit2 <- eBayes(fit[keep,], trend=TRUE) 
# hist(fit2$Amean)
# plotSA(fit2)

contrasts <- makeContrasts(FvsM = groupfemale - groupmale, 
                           levels=colnames(mm))
# Running contrast analysis
vfit <- contrasts.fit(fit, contrasts = contrasts)
# Looking at N of DEGs with adj. p <0.01 and log2FC>2
sumTable <- summary(decideTests(vfit, adjust.method = "BH", p.value = 0.05, lfc = 1))
sumTable
# Computing differential expression based on the empirical Bayes moderation of
# the standard errors towards a common value. Robust = should the estimation of
# the empirical Bayes prior parameters be robustified against outlier sample
# variances?
veBayesFit <- eBayes(vfit, robust=TRUE)
png(filename ="FIGURES/voomFits/deciduas/mat_Age_Gravity_Lane_PC1_PC2/decidua_voomFit_geneLevel_FinalModel.png",width=6,height=6,units="in",res=1200)
plotSA(veBayesFit, main = "Final model: Mean-variance trend")
dev.off()
dev.off()

# Getting the DEGs. Log2FC of 1 is equivalent to linear fold change of 2.
# Getting summary statistics for all genes
allComparisons <- colnames(contrasts)
coef = 1
#df <- data.frame(veBayesFit$sigma, veBayesFit$genes)
#df$sqrtSigma <- sqrt(df$veBayesFit.sigma)

for (i in allComparisons){
  vTopTableAll <- topTable(veBayesFit, coef=coef, n=Inf, p.value=1, lfc=0)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2/DEGs_deciduas_genes", i, "_fdr1_lfc0.txt", sep = "")
  write.table(vTopTableAll, path, sep = "\t")
  
  # Adj.p<0.05, log2FC>0
  vTopTableFDR.05 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.05, lfc=0)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2/DEGs_deciduas_genes", i, "_fdr05_lfc0.txt", sep = "")
  write.table(vTopTableFDR.05, path, sep = "\t")
  
  # log2FC>1
  vTopTableFC1 <- topTable(veBayesFit, coef=coef, n=Inf, lfc=1)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2/DEGs_deciduas_genes", i, "_lfc1.txt", sep = "")
  write.table(vTopTableFC1, path, sep = "\t")
  
  # log2FC>2
  vTopTableFC2 <- topTable(veBayesFit, coef=coef, n=Inf, lfc=2)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2/DEGs_deciduas_genes", i, "_lfc2.txt", sep = "")
  write.table(vTopTableFC2, path, sep = "\t")
  
  # Adj.p<0.05, log2FC>1
  vTopTableFDR.05FC1 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.05, lfc=1)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2/DEGs_deciduas_genes", i, "_fdr05_lfc1.txt", sep = "")
  write.table(vTopTableFDR.05FC1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>0
  vTopTableFDR.01 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=0)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2/DEGs_deciduas_genes", i, "_fdr001_lfc0.txt", sep = "")
  write.table(vTopTableFDR.01, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>1
  vTopTableFDR.01FC1 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=1)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2/DEGs_deciduas_genes", i, "_fdr001_lfc1.txt", sep = "")
  write.table(vTopTableFDR.01FC1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>2
  vTopTableFDR.01FC2 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=2)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2/DEGs_deciduas_genes", i, "_fdr001_lfc2.txt", sep = "")
  write.table(vTopTableFDR.01FC2, path, sep = "\t")
  
  coef = coef+1
}

top.table <- topTable(veBayesFit, sort.by = "P", n = Inf)
head(top.table, 10)
length(which(top.table$adj.P.Val < 0.05))
summary(decideTests(veBayesFit))
dt <- decideTests(veBayesFit)

MvsF_All <- topTable(veBayesFit, n=Inf, coef=1)
P.Value <- c(MvsF_All$P.Value)
logFC <- c(MvsF_All$logFC)
adj.P.Val <- c(MvsF_All$adj.P.Val)
df <- data.frame(P.Value, logFC, adj.P.Val, MvsF_All$Geneid, MvsF_All$chr)

# Male bias
df.MB <- subset(df, adj.P.Val < 0.05 & logFC < -1 & MvsF_All$chr != "chrY") #define male-biased, black
df.MB <- cbind(df.MB, rep(3, nrow(df.MB)))
colnames(df.MB)[6] <- "Color"

# not signigicant 
df.N <- subset(df, adj.P.Val >= 0.05) #define non-significant, grey
df.N <- cbind(df.N, rep(1, nrow(df.N)))
colnames(df.N)[6] <- "Color"

df.S <- subset(df, (adj.P.Val < 0.05 & logFC > -1) | (adj.P.Val < 0.05 & logFC < 1)) #define non-significant, grey
df.S <- cbind(df.S, rep(2, nrow(df.S)))
colnames(df.S)[6] <- "Color"

# Female bias
df.FB <- subset(df,  adj.P.Val < 0.05 & logFC > 1 & MvsF_All$chr != "chrX") #define female-biased, black
df.FB <- cbind(df.FB, rep(3, nrow(df.FB)))
colnames(df.FB)[6] <- "Color"

# X-linked and significant 
df.X <- subset(df, (adj.P.Val < 0.05 & MvsF_All$chr == "chrX" & logFC > 1) | (adj.P.Val < 0.05 & MvsF_All$chr == "chrX" & logFC < -1)) #define X-chromosomal, red
df.X <- cbind(df.X, rep(4, nrow(df.X)))
colnames(df.X)[6] <- "Color"

# Y-linked and signigicant 
df.Y <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr == "chrY" & logFC < -1) #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(5, nrow(df.Y)))
colnames(df.Y)[6] <- "Color"

df.t <- rbind(df.N, df.S, df.MB, df.FB, df.X, df.Y)
df.t$Color <- as.factor(df.t$Color)
colnames(df.t)[4] <- "Geneid"
colnames(df.t)[5] <- "chr"
colnames(df.t)[6] <- "color"
dim(df.t)

dim(df)
dim(df.N)
dim(df.S)
dim(df.MB)
dim(df.Y)
dim(df.FB)
dim(df.X)

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 65)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") 
                      | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") 
                      | adj.P.Val<0.05 & logFC < -1 & (chr=="chrX") & (color=="4") 
                      |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))

#subset_data <- subset(df.t, (logFC < -1.5 & P.Value <0.01 & color=="1") | (logFC > 1.5 & P.Value <0.01 & color=="1"))


# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2/deciduas_geneLevel_fdr05_lfc1_adjpsizeSet.png", width= 480, height= 580)
p + ggtitle("decidua differential expression gene level")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  #  xlim(c(-10, 10)) + ylim(c(0, 40)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") 
                      | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") 
                      | adj.P.Val<0.05 & logFC < -1 & (chr=="chrX") & (color=="4") 
                      |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2/deciduas_geneLevel_fdr05_lfc1.png", width= 480, height= 580)
p + ggtitle("decidua differential expression gene level")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

MvsF_All <- topTable(veBayesFit, n=Inf, coef=1)
P.Value <- c(MvsF_All$P.Value)
logFC <- c(MvsF_All$logFC)
adj.P.Val <- c(MvsF_All$adj.P.Val)
df <- data.frame(P.Value, logFC, adj.P.Val, MvsF_All$Geneid, MvsF_All$chr)

# Male bias
df.MB <- subset(df, logFC < -1 & MvsF_All$chr != "chrY") #define male-biased, black
df.MB <- cbind(df.MB, rep(3, nrow(df.MB)))
colnames(df.MB)[6] <- "Color"

# not signigicant 
df.N <- subset(df, adj.P.Val >= 0.05) #define non-significant, grey
df.N <- cbind(df.N, rep(1, nrow(df.N)))
colnames(df.N)[6] <- "Color"

df.S <- subset(df, logFC > -1 | logFC < 1) #define non-significant, grey
df.S <- cbind(df.S, rep(2, nrow(df.S)))
colnames(df.S)[6] <- "Color"

# Female bias
df.FB <- subset(df, logFC > 1 & MvsF_All$chr != "chrX") #define female-biased, black
df.FB <- cbind(df.FB, rep(3, nrow(df.FB)))
colnames(df.FB)[6] <- "Color"

# X-linked and significant 
df.X <- subset(df, (MvsF_All$chr == "chrX" & logFC > 1) | (MvsF_All$chr == "chrX" & logFC < -1)) #define X-chromosomal, red
df.X <- cbind(df.X, rep(4, nrow(df.X)))
colnames(df.X)[6] <- "Color"

# Y-linked and signigicant 
df.Y <- subset(df,  MvsF_All$chr == "chrY" & logFC < -1) #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(5, nrow(df.Y)))
colnames(df.Y)[6] <- "Color"

df.t <- rbind(df.N, df.S, df.MB, df.FB, df.X, df.Y)
df.t$Color <- as.factor(df.t$Color)
colnames(df.t)[4] <- "Geneid"
colnames(df.t)[5] <- "chr"
colnames(df.t)[6] <- "color"
dim(df.t)

dim(df)
dim(df.N)
dim(df.S)
dim(df.MB)
dim(df.Y)
dim(df.FB)
dim(df.X)

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 65)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, logFC < -1 & (chr=="chrY") & (color=="5") 
                      |logFC > 1 & (chr=="chrX") & (color=="4") 
                      | logFC < -1 & (chr=="chrX") & (color=="4") 
                      | (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2/deciduas_geneLevel_lfc1_sizeSet.png", width= 480, height= 580)
p + ggtitle("decidua differential expression gene level")+
  #geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  #  xlim(c(-10, 10)) + ylim(c(0, 40)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, logFC < -1 & (chr=="chrY") & (color=="5") 
                      | logFC > 1 & (chr=="chrX") & (color=="4") 
                      | logFC < -1 & (chr=="chrX") & (color=="4") 
                      | (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2/deciduas_geneLevel_lfc1.png", width= 480, height= 580)
p + ggtitle("decidua differential expression gene level")+
  #geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

MvsF_All <- topTable(veBayesFit, n=Inf, coef=1)
P.Value <- c(MvsF_All$P.Value)
logFC <- c(MvsF_All$logFC)
adj.P.Val <- c(MvsF_All$adj.P.Val)
df <- data.frame(P.Value, logFC, adj.P.Val, MvsF_All$Geneid, MvsF_All$chr)

# Male bias
df.MB <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr != "chrY") #define male-biased, black
df.MB <- cbind(df.MB, rep(3, nrow(df.MB)))
colnames(df.MB)[6] <- "Color"

# not signigicant 
df.N <- subset(df, adj.P.Val >= 0.05) #define non-significant, grey
df.N <- cbind(df.N, rep(1, nrow(df.N)))
colnames(df.N)[6] <- "Color"

df.S <- subset(df, adj.P.Val < 0.05) #define non-significant, grey
df.S <- cbind(df.S, rep(2, nrow(df.S)))
colnames(df.S)[6] <- "Color"

# Female bias
df.FB <- subset(df,  adj.P.Val < 0.05 & MvsF_All$chr != "chrX") #define female-biased, black
df.FB <- cbind(df.FB, rep(3, nrow(df.FB)))
colnames(df.FB)[6] <- "Color"

# X-linked and significant 
df.X <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr == "chrX" ) #define X-chromosomal, red
df.X <- cbind(df.X, rep(4, nrow(df.X)))
colnames(df.X)[6] <- "Color"

# Y-linked and signigicant 
df.Y <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr == "chrY") #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(5, nrow(df.Y)))
colnames(df.Y)[6] <- "Color"

df.t <- rbind(df.N, df.S, df.MB, df.FB, df.X, df.Y)
df.t$Color <- as.factor(df.t$Color)
colnames(df.t)[4] <- "Geneid"
colnames(df.t)[5] <- "chr"
colnames(df.t)[6] <- "color"
dim(df.t)

dim(df)
dim(df.N)
dim(df.S)
dim(df.MB)
dim(df.Y)
dim(df.FB)
dim(df.X)

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 65)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, adj.P.Val<0.05 & (chr=="chrY") & (color=="5") 
                      | adj.P.Val<0.05 & (chr=="chrX") & (color=="4") 
                      | adj.P.Val<0.05 & (chr=="chrX") & (color=="4") 
                      |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2/deciduas_geneLevel_fdr05_sizeSet.png", width= 480, height= 580)
p + ggtitle("decidua differential expression gene level")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  #  xlim(c(-10, 10)) + ylim(c(0, 40)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, adj.P.Val<0.05  & (chr=="chrY") & (color=="5") 
                      | adj.P.Val<0.05  & (chr=="chrX") & (color=="4") 
                      | adj.P.Val<0.05  & (chr=="chrX") & (color=="4") 
                      |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2/deciduas_geneLevel_fdr05.png", width= 480, height= 580)
p + ggtitle("decidua differential expression gene level")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

#------------- VOOM 
#--------- tran level Voom 
group <-interaction(DGE_tran$samples$Offspring_sex) #, DGE$samples$batch, DGE$samples$Lane, DGE$samples$Reported_race) #, DGE$samples$batch, DGE$samples$Reported_race, DGE$samples$site, DGE$samples$Lane, DGE$samples)
mm <- model.matrix(~0 + group + Maternal_age + Gravidity + Lane + PC1 + PC2)
#mm <- model.matrix(~0 + group + Maternal_age + Gravidity + MethodConception + Gestational_age + Lane + PC1 + PC2)
colnames(mm) <- gsub("v\\$targets\\$Offspring_sex", "", colnames(mm))
dim(mm)

y <- voom(DGE_tran, mm, plot = T) # plot T for plotting 

# Fitting linear models in limma
fit <- lmFit(y, mm)

png(filename ="FIGURES/voomFits/deciduas/mat_Age_Gravity_Lane_PC1_PC2/decidua_voomFit_trans.png",width=6,height=6,units="in",res=1200)
v <- voomWithQualityWeights(DGE_tran, mm, plot=TRUE)
dev.off()
dev.off()

voomCounts <- v$E
voomCountsMatrix <- data.matrix(voomCounts, rownames.force = NA)
#write.table(voomCounts, "decidua_voom_cpm_genelevel.txt", sep="\t", row.names = FALSE, quote = FALSE)
genes <- v$genes
genesVoom <- cbind(genes, voomCounts)

# #-------------------------------
# variance Part
# Fit model and extract results. A linear mixed model is used
# each entry in results is a regression model fit on a single gene
# form <- ~ Maternal_age + Gravidity + Parity + Pre.Pregnancy_BMI + BirthWeight + (1|MethodConception) + (1|Reported_race) + (1|Lane) + (1|Offspring_sex)
# 
# # 2) extract variance fractions from each model fit for each gene, returns fraction of variation attributable to each variable
# # Interpretation: the variance explained by each variables after correcting for all other variables
# 
# cl <- makeCluster(4)
# registerDoParallel(cl)
# 
# # Fitting the linear model
# varPart <- fitExtractVarPartModel(voomCounts, form, pheno_ExRemovals)
# 
# #varPart <- fitExtractVarPartModel(df, form, pheno_ExRemovals )
# #write.table(varPart, "variancePart/varPart_deciduas.txt")
# # to try and run this
# #the brglm2 package has a method for detecting complete separation:
# # install.packages("devtools")
# # devtools::install_github("ikosmidis/brglm2")
# # library("brglm2")
# #glm(form, data = genesVoom, pheno_ExRemovals, family = binomial, method="detect_separation")
# # data frame of pheno and counts....
# #This should say whether complete separation occurs, and in which (combinations of) variables, e.g.
# # sort variables (i.e. columns) by median fraction of variance explained
# genesVP <- cbind(genes, varPart)
# write.table(genesVP, "variancePart/deciduas/varPartgenes_deciduas_trans.txt")
# 
# # violin plot of contribution of each variable to total variance
# vp <- sortCols( varPart )
# png(filename ="variancePart/deciduas/varPart_deciduas_trans.png",width=6,height=6,units="in",res=1200)
# plotVarPart( vp )
# dev.off()
# dev.off()

#---- end var part
#-------------------------------

#write.table(genes, "decidua_genes.txt", row.names = FALSE, quote = FALSE)
#v$genes
fit <- lmFit(v, mm)
png(filename ="FIGURES/voomFits/deciduas/mat_Age_Gravity_Lane_PC1_PC2/decidua_fitAmeanHist_fpkm1_trans.png",width=6,height=6,units="in",res=1200)
hist(fit$Amean) 
dev.off()
dev.off()
plotSA(fit)
# keep <- fit$Amean > 0
# fit2 <- eBayes(fit[keep,], trend=TRUE) 
# hist(fit2$Amean)
# plotSA(fit2)

contrasts <- makeContrasts(FvsM = groupfemale - groupmale, 
                           levels=colnames(mm))
# Running contrast analysis
vfit <- contrasts.fit(fit, contrasts = contrasts)
# Looking at N of DEGs with adj. p <0.01 and log2FC>2
sumTable <- summary(decideTests(vfit, adjust.method = "BH", p.value = 0.05, lfc = 1))
sumTable
# Computing differential expression based on the empirical Bayes moderation of
# the standard errors towards a common value. Robust = should the estimation of
# the empirical Bayes prior parameters be robustified against outlier sample
# variances?
veBayesFit <- eBayes(vfit, robust=TRUE)
png(filename ="FIGURES/voomFits/deciduas/mat_Age_Gravity_Lane_PC1_PC2/decidua_voomFit_geneLevel_FinalModel_trans.png",width=6,height=6,units="in",res=1200)
plotSA(veBayesFit, main = "Final model: Mean-variance trend")
dev.off()
dev.off()

# Getting the DEGs. Log2FC of 1 is equivalent to linear fold change of 2.
# Getting summary statistics for all genes
allComparisons <- colnames(contrasts)
coef = 1

for (i in allComparisons){
  vTopTableAll <- topTable(veBayesFit, coef=coef, n=Inf, p.value=1, lfc=0)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2/DEGs_deciduas_trans", i, "_fdr1_lfc0.txt", sep = "")
  write.table(vTopTableAll, path, sep = "\t")
  
  # Adj.p<0.05, log2FC>0
  vTopTableFDR.05 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.05, lfc=0)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2/DEGs_deciduas_trans", i, "_fdr05_lfc0.txt", sep = "")
  write.table(vTopTableFDR.05, path, sep = "\t")
  
  # log2FC>1
  vTopTableFC1 <- topTable(veBayesFit, coef=coef, n=Inf, lfc=1)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2/DEGs_deciduas_trans", i, "_lfc1.txt", sep = "")
  write.table(vTopTableFC1, path, sep = "\t")
  
  # log2FC>2
  vTopTableFC2 <- topTable(veBayesFit, coef=coef, n=Inf, lfc=2)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2/DEGs_deciduas_trans", i, "_lfc2.txt", sep = "")
  write.table(vTopTableFC2, path, sep = "\t")
  
  # Adj.p<0.05, log2FC>1
  vTopTableFDR.05FC1 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.05, lfc=1)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2/DEGs_deciduas_trans", i, "_fdr05_lfc1.txt", sep = "")
  write.table(vTopTableFDR.05FC1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>0
  vTopTableFDR.01 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=0)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2/DEGs_deciduas_trans", i, "_fdr001_lfc0.txt", sep = "")
  write.table(vTopTableFDR.01, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>1
  vTopTableFDR.01FC1 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=1)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2/DEGs_deciduas_trans", i, "_fdr001_lfc1.txt", sep = "")
  write.table(vTopTableFDR.01FC1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>2
  vTopTableFDR.01FC2 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=2)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2/DEGs_deciduas_trans", i, "_fdr001_lfc2.txt", sep = "")
  write.table(vTopTableFDR.01FC2, path, sep = "\t")
  
  coef = coef+1
}

top.table <- topTable(veBayesFit, sort.by = "P", n = Inf)
head(top.table, 10)
length(which(top.table$adj.P.Val < 0.05))
summary(decideTests(veBayesFit))
dt <- decideTests(veBayesFit)

MvsF_All <- topTable(veBayesFit, n=Inf, coef=1)
P.Value <- c(MvsF_All$P.Value)
logFC <- c(MvsF_All$logFC)
adj.P.Val <- c(MvsF_All$adj.P.Val)
df <- data.frame(P.Value, logFC, adj.P.Val, MvsF_All$Geneid, MvsF_All$chr, MvsF_All$gene_name)

# Male bias
df.MB <- subset(df, adj.P.Val < 0.05 & logFC < -1 & MvsF_All$chr != "chrY") #define male-biased, black
df.MB <- cbind(df.MB, rep(3, nrow(df.MB)))
colnames(df.MB)[7] <- "color"

# not signigicant 
df.N <- subset(df, adj.P.Val >= 0.05) #define non-significant, grey
df.N <- cbind(df.N, rep(1, nrow(df.N)))
colnames(df.N)[7] <- "color"

df.S <- subset(df, (adj.P.Val < 0.05 & logFC > -1) | (adj.P.Val < 0.05 & logFC < 1)) #define non-significant, grey
df.S <- cbind(df.S, rep(2, nrow(df.S)))
colnames(df.S)[7] <- "color"

# Female bias
df.FB <- subset(df,  adj.P.Val < 0.05 & logFC > 1 & MvsF_All$chr != "chrX") #define female-biased, black
df.FB <- cbind(df.FB, rep(3, nrow(df.FB)))
colnames(df.FB)[7] <- "color"

# X-linked and significant 
df.X <- subset(df, (adj.P.Val < 0.05 & MvsF_All$chr == "chrX" & logFC > 1) | (adj.P.Val < 0.05 & MvsF_All$chr == "chrX" & logFC < -1)) #define X-chromosomal, red
df.X <- cbind(df.X, rep(4, nrow(df.X)))
colnames(df.X)[7] <- "color"

# Y-linked and signigicant 
df.Y <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr == "chrY" & logFC < -1) #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(5, nrow(df.Y)))
colnames(df.Y)[7] <- "color"

df.t <- rbind(df.N, df.S, df.MB, df.FB, df.X, df.Y)
df.t$color <- as.factor(df.t$color)
colnames(df.t)[4] <- "Geneid"
colnames(df.t)[5] <- "chr"
colnames(df.t)[6] <- "gene_name"
colnames(df.t)[7] <- "color"
dim(df.t)

dim(df)
dim(df.N)
dim(df.S)
dim(df.MB)
dim(df.Y)
dim(df.FB)
dim(df.X)

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 65)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") 
                      | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") 
                      | adj.P.Val<0.05 & logFC < -1 & (chr=="chrX") & (color=="4") 
                      |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))

#subset_data <- subset(df.t, (logFC < -1.5 & P.Value <0.01 & color=="1") | (logFC > 1.5 & P.Value <0.01 & color=="1"))


# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2/deciduas_tranLevel_fdr05_lfc1_adjpsizeSet.png", width= 480, height= 580)
p + ggtitle("decidua differential expression transcript level")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene_name), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  #  xlim(c(-10, 10)) + ylim(c(0, 40)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") 
                      | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") 
                      | adj.P.Val<0.05 & logFC < -1 & (chr=="chrX") & (color=="4") 
                      |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2/deciduas_tranLevel_fdr05_lfc1.png", width= 480, height= 580)
p + ggtitle("decidua differential expression transcript level")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene_name), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

MvsF_All <- topTable(veBayesFit, n=Inf, coef=1)
P.Value <- c(MvsF_All$P.Value)
logFC <- c(MvsF_All$logFC)
adj.P.Val <- c(MvsF_All$adj.P.Val)
df <- data.frame(P.Value, logFC, adj.P.Val, MvsF_All$Geneid, MvsF_All$chr, MvsF_All$gene_name)

# Male bias
df.MB <- subset(df, logFC < -1 & MvsF_All$chr != "chrY") #define male-biased, black
df.MB <- cbind(df.MB, rep(3, nrow(df.MB)))
colnames(df.MB)[7] <- "color"

# not signigicant 
df.N <- subset(df, adj.P.Val >= 0.05) #define non-significant, grey
df.N <- cbind(df.N, rep(1, nrow(df.N)))
colnames(df.N)[7] <- "color"

df.S <- subset(df, logFC > -1 | logFC < 1) #define non-significant, grey
df.S <- cbind(df.S, rep(2, nrow(df.S)))
colnames(df.S)[7] <- "color"

# Female bias
df.FB <- subset(df, logFC > 1 & MvsF_All$chr != "chrX") #define female-biased, black
df.FB <- cbind(df.FB, rep(3, nrow(df.FB)))
colnames(df.FB)[7] <- "color"

# X-linked and significant 
df.X <- subset(df, (MvsF_All$chr == "chrX" & logFC > 1) | (MvsF_All$chr == "chrX" & logFC < -1)) #define X-chromosomal, red
df.X <- cbind(df.X, rep(4, nrow(df.X)))
colnames(df.X)[7] <- "color"

# Y-linked and signigicant 
df.Y <- subset(df,  MvsF_All$chr == "chrY" & logFC < -1) #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(5, nrow(df.Y)))
colnames(df.Y)[7] <- "color"

df.t <- rbind(df.N, df.S, df.MB, df.FB, df.X, df.Y)
df.t$color <- as.factor(df.t$color)
colnames(df.t)[4] <- "Geneid"
colnames(df.t)[5] <- "chr"
colnames(df.t)[6] <- "gene_name"
colnames(df.t)[7] <- "color"
dim(df.t)

dim(df)
dim(df.N)
dim(df.S)
dim(df.MB)
dim(df.Y)
dim(df.FB)
dim(df.X)

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 65)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, logFC < -1 & (chr=="chrY") & (color=="5") 
                      |logFC > 1 & (chr=="chrX") & (color=="4") 
                      | logFC < -1 & (chr=="chrX") & (color=="4") 
                      | (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2/deciduas_tranLevel_lfc1_sizeSet.png", width= 480, height= 580)
p + ggtitle("decidua differential expression transcript level")+
  # geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene_name), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  #  xlim(c(-10, 10)) + ylim(c(0, 40)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, logFC < -1 & (chr=="chrY") & (color=="5") 
                      | logFC > 1 & (chr=="chrX") & (color=="4") 
                      | logFC < -1 & (chr=="chrX") & (color=="4") 
                      | (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2/deciduas_tranLevel_lfc1.png", width= 480, height= 580)
p + ggtitle("decidua differential expression transcript level")+
  # geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene_name), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

MvsF_All <- topTable(veBayesFit, n=Inf, coef=1)
P.Value <- c(MvsF_All$P.Value)
logFC <- c(MvsF_All$logFC)
adj.P.Val <- c(MvsF_All$adj.P.Val)
df <- data.frame(P.Value, logFC, adj.P.Val, MvsF_All$Geneid, MvsF_All$chr, MvsF_All$gene_name)

# Male bias
df.MB <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr != "chrY") #define male-biased, black
df.MB <- cbind(df.MB, rep(3, nrow(df.MB)))
colnames(df.MB)[7] <- "color"

# not signigicant 
df.N <- subset(df, adj.P.Val >= 0.05) #define non-significant, grey
df.N <- cbind(df.N, rep(1, nrow(df.N)))
colnames(df.N)[7] <- "color"

df.S <- subset(df, adj.P.Val < 0.05) #define non-significant, grey
df.S <- cbind(df.S, rep(2, nrow(df.S)))
colnames(df.S)[7] <- "color"

# Female bias
df.FB <- subset(df,  adj.P.Val < 0.05 & MvsF_All$chr != "chrX") #define female-biased, black
df.FB <- cbind(df.FB, rep(3, nrow(df.FB)))
colnames(df.FB)[7] <- "color"

# X-linked and significant 
df.X <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr == "chrX" ) #define X-chromosomal, red
df.X <- cbind(df.X, rep(4, nrow(df.X)))
colnames(df.X)[7] <- "color"

# Y-linked and signigicant 
df.Y <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr == "chrY") #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(5, nrow(df.Y)))
colnames(df.Y)[7] <- "color"

df.t <- rbind(df.N, df.S, df.MB, df.FB, df.X, df.Y)
df.t$color <- as.factor(df.t$color)
colnames(df.t)[4] <- "Geneid"
colnames(df.t)[5] <- "chr"
colnames(df.t)[6] <- "gene_name"
colnames(df.t)[7] <- "color"
dim(df.t)

dim(df)
dim(df.N)
dim(df.S)
dim(df.MB)
dim(df.Y)
dim(df.FB)
dim(df.X)

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 65)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, adj.P.Val<0.05 & (chr=="chrY") & (color=="5") 
                      | adj.P.Val<0.05 & (chr=="chrX") & (color=="4") 
                      | adj.P.Val<0.05 & (chr=="chrX") & (color=="4") 
                      |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2/deciduas_tranLevel_fdr05_sizeSet.png", width= 480, height= 580)
p + ggtitle("decidua differential expression transcript level")+
  #  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene_name), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  #  xlim(c(-10, 10)) + ylim(c(0, 40)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, adj.P.Val<0.05  & (chr=="chrY") & (color=="5") 
                      | adj.P.Val<0.05  & (chr=="chrX") & (color=="4") 
                      | adj.P.Val<0.05  & (chr=="chrX") & (color=="4") 
                      |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2/deciduas_tranLevel_fdr05.png", width= 480, height= 580)
p + ggtitle("decidua differential expression transcript level")+
  #geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene_name), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------

#------------- VOOM 
#--------- gene level Voom 
group <-interaction(DGE_gene$samples$Offspring_sex) #, DGE$samples$batch, DGE$samples$Lane, DGE$samples$Reported_race) #, DGE$samples$batch, DGE$samples$Reported_race, DGE$samples$site, DGE$samples$Lane, DGE$samples)
mm <- model.matrix(~0 + group + Maternal_age + Gravidity + Gestational_age + Lane + PC1 + PC2)
colnames(mm) <- gsub("v\\$targets\\$Offspring_sex", "", colnames(mm))
dim(mm)

y <- voom(DGE_gene, mm, plot = T) # plot T for plotting 

# Fitting linear models in limma
fit <- lmFit(y, mm)

png(filename ="FIGURES/voomFits/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/decidua_voomFit.png",width=6,height=6,units="in",res=1200)
v <- voomWithQualityWeights(DGE_gene, mm, plot=TRUE)
dev.off()
dev.off()

voomCounts <- v$E
voomCountsMatrix <- data.matrix(voomCounts, rownames.force = NA)
#write.table(voomCounts, "decidua_voom_cpm_genelevel.txt", sep="\t", row.names = FALSE, quote = FALSE)
genes <- v$genes
genesVoom <- cbind(genes, voomCounts)


#write.table(genes, "decidua_genes.txt", row.names = FALSE, quote = FALSE)
#v$genes
fit <- lmFit(v, mm)
png(filename ="FIGURES/voomFits/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/decidua_fitAmeanHist_fpkm1.png",width=6,height=6,units="in",res=1200)
hist(fit$Amean) 
dev.off()
dev.off()
plotSA(fit)
# keep <- fit$Amean > 0
# fit2 <- eBayes(fit[keep,], trend=TRUE) 
# hist(fit2$Amean)
# plotSA(fit2)

contrasts <- makeContrasts(FvsM = groupfemale - groupmale, 
                           levels=colnames(mm))
# Running contrast analysis
vfit <- contrasts.fit(fit, contrasts = contrasts)
# Looking at N of DEGs with adj. p <0.01 and log2FC>2
sumTable <- summary(decideTests(vfit, adjust.method = "BH", p.value = 0.05, lfc = 1))
sumTable

may <- cbind(genes, voomCounts)
may$length <- NULL
may2 <- melt(may, id=c("Geneid", "chr"))
df_merged <- merge(may2, pheno_ExRemovals, by.x="variable", by.y="sample", all.x=TRUE)
write.table(df_merged, "genelists/deciduas/vfit_decidua_df_genes.txt", sep = "\t")


# Computing differential expression based on the empirical Bayes moderation of
# the standard errors towards a common value. Robust = should the estimation of
# the empirical Bayes prior parameters be robustified against outlier sample
# variances?
veBayesFit <- eBayes(vfit, robust=TRUE)
png(filename ="FIGURES/voomFits/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/decidua_voomFit_geneLevel_FinalModel.png",width=6,height=6,units="in",res=1200)
plotSA(veBayesFit, main = "Final model: Mean-variance trend")
dev.off()
dev.off()

# Getting the DEGs. Log2FC of 1 is equivalent to linear fold change of 2.
# Getting summary statistics for all genes
allComparisons <- colnames(contrasts)
coef = 1
#df <- data.frame(veBayesFit$sigma, veBayesFit$genes)
#df$sqrtSigma <- sqrt(df$veBayesFit.sigma)

for (i in allComparisons){
  vTopTableAll <- topTable(veBayesFit, coef=coef, n=Inf, p.value=1, lfc=0)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/DEGs_deciduas_genes", i, "_fdr1_lfc0.txt", sep = "")
  write.table(vTopTableAll, path, sep = "\t")
  
  # Adj.p<0.05, log2FC>0
  vTopTableFDR.05 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.05, lfc=0)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/DEGs_deciduas_genes", i, "_fdr05_lfc0.txt", sep = "")
  write.table(vTopTableFDR.05, path, sep = "\t")
  
  # log2FC>1
  vTopTableFC1 <- topTable(veBayesFit, coef=coef, n=Inf, lfc=1)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/DEGs_deciduas_genes", i, "_lfc1.txt", sep = "")
  write.table(vTopTableFC1, path, sep = "\t")
  
  # log2FC>2
  vTopTableFC2 <- topTable(veBayesFit, coef=coef, n=Inf, lfc=2)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/DEGs_deciduas_genes", i, "_lfc2.txt", sep = "")
  write.table(vTopTableFC2, path, sep = "\t")
  
  # Adj.p<0.05, log2FC>1
  vTopTableFDR.05FC1 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.05, lfc=1)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/DEGs_deciduas_genes", i, "_fdr05_lfc1.txt", sep = "")
  write.table(vTopTableFDR.05FC1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>0
  vTopTableFDR.01 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=0)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/DEGs_deciduas_genes", i, "_fdr001_lfc0.txt", sep = "")
  write.table(vTopTableFDR.01, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>1
  vTopTableFDR.01FC1 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=1)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/DEGs_deciduas_genes", i, "_fdr001_lfc1.txt", sep = "")
  write.table(vTopTableFDR.01FC1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>2
  vTopTableFDR.01FC2 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=2)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/DEGs_deciduas_genes", i, "_fdr001_lfc2.txt", sep = "")
  write.table(vTopTableFDR.01FC2, path, sep = "\t")
  
  coef = coef+1
}

top.table <- topTable(veBayesFit, sort.by = "P", n = Inf)
head(top.table, 10)
length(which(top.table$adj.P.Val < 0.05))
summary(decideTests(veBayesFit))
dt <- decideTests(veBayesFit)

MvsF_All <- topTable(veBayesFit, n=Inf, coef=1)
P.Value <- c(MvsF_All$P.Value)
logFC <- c(MvsF_All$logFC)
adj.P.Val <- c(MvsF_All$adj.P.Val)
df <- data.frame(P.Value, logFC, adj.P.Val, MvsF_All$Geneid, MvsF_All$chr)

# Male bias
df.MB <- subset(df, adj.P.Val < 0.05 & logFC < -1 & MvsF_All$chr != "chrY") #define male-biased, black
df.MB <- cbind(df.MB, rep(3, nrow(df.MB)))
colnames(df.MB)[6] <- "Color"

# not signigicant 
df.N <- subset(df, adj.P.Val >= 0.05) #define non-significant, grey
df.N <- cbind(df.N, rep(1, nrow(df.N)))
colnames(df.N)[6] <- "Color"

df.S <- subset(df, (adj.P.Val < 0.05 & logFC > -1) | (adj.P.Val < 0.05 & logFC < 1)) #define non-significant, grey
df.S <- cbind(df.S, rep(2, nrow(df.S)))
colnames(df.S)[6] <- "Color"

# Female bias
df.FB <- subset(df,  adj.P.Val < 0.05 & logFC > 1 & MvsF_All$chr != "chrX") #define female-biased, black
df.FB <- cbind(df.FB, rep(3, nrow(df.FB)))
colnames(df.FB)[6] <- "Color"

# X-linked and significant 
df.X <- subset(df, (adj.P.Val < 0.05 & MvsF_All$chr == "chrX" & logFC > 1) | (adj.P.Val < 0.05 & MvsF_All$chr == "chrX" & logFC < -1)) #define X-chromosomal, red
df.X <- cbind(df.X, rep(4, nrow(df.X)))
colnames(df.X)[6] <- "Color"

# Y-linked and signigicant 
df.Y <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr == "chrY" & logFC < -1) #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(5, nrow(df.Y)))
colnames(df.Y)[6] <- "Color"

df.t <- rbind(df.N, df.S, df.MB, df.FB, df.X, df.Y)
df.t$Color <- as.factor(df.t$Color)
colnames(df.t)[4] <- "Geneid"
colnames(df.t)[5] <- "chr"
colnames(df.t)[6] <- "color"
dim(df.t)

dim(df)
dim(df.N)
dim(df.S)
dim(df.MB)
dim(df.Y)
dim(df.FB)
dim(df.X)

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 65)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") 
                      | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") 
                      | adj.P.Val<0.05 & logFC < -1 & (chr=="chrX") & (color=="4") 
                      |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))

#subset_data <- subset(df.t, (logFC < -1.5 & P.Value <0.01 & color=="1") | (logFC > 1.5 & P.Value <0.01 & color=="1"))


# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/deciduas_geneLevel_fdr05_lfc1_adjpsizeSet.png", width= 480, height= 580)
p + ggtitle("decidua differential expression gene level")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  xlim(c(-7, 7)) +
  #ylim(c(0, 40)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=18), 
        axis.text.x=element_text(size=18)) +
  theme(axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=18))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") 
                      | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") 
                      | adj.P.Val<0.05 & logFC < -1 & (chr=="chrX") & (color=="4") 
                      |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))

subset_data <- subset(df.t, logFC < -1 & (chr=="chrY") & (color=="5") 
                      | logFC > 1 & (chr=="chrX") & (color=="4") 
                      | logFC < -1 & (chr=="chrX") & (color=="4") 
                      | (chr!="chrX") & (chr!="chrY") & (color=="3") & logFC >1
                      | (chr!="chrX") & (chr!="chrY") & (color=="3") & logFC < -1 )
subset_data<- unique(subset_data)

# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/deciduas_geneLevel_fdr05_lfc1.png", width= 480, height= 580)
p2 <- p + ggtitle("")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), 
                  segment.alpha = 1, size = 5.5, fontface="italic") + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
p2
dev.off()
dev.off()


MvsF_All <- topTable(veBayesFit, n=Inf, coef=1)
P.Value <- c(MvsF_All$P.Value)
logFC <- c(MvsF_All$logFC)
adj.P.Val <- c(MvsF_All$adj.P.Val)
df <- data.frame(P.Value, logFC, adj.P.Val, MvsF_All$Geneid, MvsF_All$chr)

# Male bias
df.MB <- subset(df, logFC < -1 & MvsF_All$chr != "chrY") #define male-biased, black
df.MB <- cbind(df.MB, rep(3, nrow(df.MB)))
colnames(df.MB)[6] <- "Color"

# not signigicant 
df.N <- subset(df, adj.P.Val >= 0.05) #define non-significant, grey
df.N <- cbind(df.N, rep(1, nrow(df.N)))
colnames(df.N)[6] <- "Color"

df.S <- subset(df, logFC > -1 | logFC < 1) #define non-significant, grey
df.S <- cbind(df.S, rep(2, nrow(df.S)))
colnames(df.S)[6] <- "Color"

# Female bias
df.FB <- subset(df, logFC > 1 & MvsF_All$chr != "chrX") #define female-biased, black
df.FB <- cbind(df.FB, rep(3, nrow(df.FB)))
colnames(df.FB)[6] <- "Color"

# X-linked and significant 
df.X <- subset(df, (MvsF_All$chr == "chrX" & logFC > 1) | (MvsF_All$chr == "chrX" & logFC < -1)) #define X-chromosomal, red
df.X <- cbind(df.X, rep(4, nrow(df.X)))
colnames(df.X)[6] <- "Color"

# Y-linked and signigicant 
df.Y <- subset(df,  MvsF_All$chr == "chrY" & logFC < -1) #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(5, nrow(df.Y)))
colnames(df.Y)[6] <- "Color"

df.t <- rbind(df.N, df.S, df.MB, df.FB, df.X, df.Y)
df.t$Color <- as.factor(df.t$Color)
colnames(df.t)[4] <- "Geneid"
colnames(df.t)[5] <- "chr"
colnames(df.t)[6] <- "color"
dim(df.t)

dim(df)
dim(df.N)
dim(df.S)
dim(df.MB)
dim(df.Y)
dim(df.FB)
dim(df.X)

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 65)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, logFC < -1 & (chr=="chrY") & (color=="5") 
                      |logFC > 1 & (chr=="chrX") & (color=="4") 
                      | logFC < -1 & (chr=="chrX") & (color=="4") 
                      | (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/deciduas_geneLevel_lfc1_sizeSet.png", width= 480, height= 580)
p + ggtitle("decidua differential expression gene level")+
  #geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  #  xlim(c(-10, 10)) + ylim(c(0, 40)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, logFC < -1 & (chr=="chrY") & (color=="5") 
                      | logFC > 1 & (chr=="chrX") & (color=="4") 
                      | logFC < -1 & (chr=="chrX") & (color=="4") 
                      | (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/deciduas_geneLevel_lfc1.png", width= 480, height= 580)
p + ggtitle("decidua differential expression gene level")+
  #geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

MvsF_All <- topTable(veBayesFit, n=Inf, coef=1)
P.Value <- c(MvsF_All$P.Value)
logFC <- c(MvsF_All$logFC)
adj.P.Val <- c(MvsF_All$adj.P.Val)
df <- data.frame(P.Value, logFC, adj.P.Val, MvsF_All$Geneid, MvsF_All$chr)

# Male bias
df.MB <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr != "chrY") #define male-biased, black
df.MB <- cbind(df.MB, rep(3, nrow(df.MB)))
colnames(df.MB)[6] <- "Color"

# not signigicant 
df.N <- subset(df, adj.P.Val >= 0.05) #define non-significant, grey
df.N <- cbind(df.N, rep(1, nrow(df.N)))
colnames(df.N)[6] <- "Color"

df.S <- subset(df, adj.P.Val < 0.05) #define non-significant, grey
df.S <- cbind(df.S, rep(2, nrow(df.S)))
colnames(df.S)[6] <- "Color"

# Female bias
df.FB <- subset(df,  adj.P.Val < 0.05 & MvsF_All$chr != "chrX") #define female-biased, black
df.FB <- cbind(df.FB, rep(3, nrow(df.FB)))
colnames(df.FB)[6] <- "Color"

# X-linked and significant 
df.X <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr == "chrX" ) #define X-chromosomal, red
df.X <- cbind(df.X, rep(4, nrow(df.X)))
colnames(df.X)[6] <- "Color"

# Y-linked and signigicant 
df.Y <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr == "chrY") #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(5, nrow(df.Y)))
colnames(df.Y)[6] <- "Color"

df.t <- rbind(df.N, df.S, df.MB, df.FB, df.X, df.Y)
df.t$Color <- as.factor(df.t$Color)
colnames(df.t)[4] <- "Geneid"
colnames(df.t)[5] <- "chr"
colnames(df.t)[6] <- "color"
dim(df.t)

dim(df)
dim(df.N)
dim(df.S)
dim(df.MB)
dim(df.Y)
dim(df.FB)
dim(df.X)

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  xlim(c(-7, 7)) + 
  #ylim(c(0, 65)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=18), 
        axis.text.x=element_text(size=18)) +
  theme(axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=18))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
# subset_data <- subset(df.t, adj.P.Val<0.05 & (chr=="chrY") & (color=="5")
#                       | adj.P.Val<0.05 & (chr=="chrX") & (color=="4")
#                       | adj.P.Val<0.05 & (chr=="chrX") & (color=="4")
#                       |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))
subset_data <- subset(df.t, logFC < -1 & (chr=="chrY") & (color=="5") 
                      | logFC > 1 & (chr=="chrX") & (color=="4") 
                      | logFC < -1 & (chr=="chrX") & (color=="4") 
                      | (chr!="chrX") & (chr!="chrY") & (color=="3") & logFC >1
                      | (chr!="chrX") & (chr!="chrY") & (color=="3") & logFC < -1 )
subset_data<- unique(subset_data)

# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/deciduas_geneLevel_fdr05_lfc1.png", width= 440, height= 540)
p2 <- p + ggtitle("")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), 
                  segment.alpha = 1, size = 5.5, fontface="italic") + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
p2
dev.off()
dev.off()
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/deciduas_geneLevel_fdr05_sizeSet.png", width= 480, height= 580)
p + ggtitle("decidua differential expression gene level")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  #  xlim(c(-10, 10)) + ylim(c(0, 40)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, adj.P.Val<0.05  & (chr=="chrY") & (color=="5") 
                      | adj.P.Val<0.05  & (chr=="chrX") & (color=="4") 
                      | adj.P.Val<0.05  & (chr=="chrX") & (color=="4") 
                      |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/deciduas_geneLevel_fdr05.png", width= 480, height= 580)
p + ggtitle("decidua differential expression gene level")+
  #geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

#------------- VOOM 
#--------- tran level Voom 
group <-interaction(DGE_tran$samples$Offspring_sex) #, DGE$samples$batch, DGE$samples$Lane, DGE$samples$Reported_race) #, DGE$samples$batch, DGE$samples$Reported_race, DGE$samples$site, DGE$samples$Lane, DGE$samples)
#mm <- model.matrix(~0 + group + Maternal_age + Gravidity + Lane + PC1 + PC2)
mm <- model.matrix(~0 + group + Maternal_age + Gravidity + Gestational_age + Lane + PC1 + PC2)
colnames(mm) <- gsub("v\\$targets\\$Offspring_sex", "", colnames(mm))
dim(mm)

y <- voom(DGE_tran, mm, plot = T) # plot T for plotting 

# Fitting linear models in limma
fit <- lmFit(y, mm)

png(filename ="FIGURES/voomFits/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/decidua_voomFit_trans.png",width=6,height=6,units="in",res=1200)
v <- voomWithQualityWeights(DGE_tran, mm, plot=TRUE)
dev.off()
dev.off()

voomCounts <- v$E
voomCountsMatrix <- data.matrix(voomCounts, rownames.force = NA)
#write.table(voomCounts, "decidua_voom_cpm_genelevel.txt", sep="\t", row.names = FALSE, quote = FALSE)
genes <- v$genes
genesVoom <- cbind(genes, voomCounts)

#write.table(genes, "decidua_genes.txt", row.names = FALSE, quote = FALSE)
#v$genes
fit <- lmFit(v, mm)
png(filename ="FIGURES/voomFits/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/decidua_fitAmeanHist_fpkm1_trans.png",width=6,height=6,units="in",res=1200)
hist(fit$Amean) 
dev.off()
dev.off()
plotSA(fit)
# keep <- fit$Amean > 0
# fit2 <- eBayes(fit[keep,], trend=TRUE) 
# hist(fit2$Amean)
# plotSA(fit2)

contrasts <- makeContrasts(FvsM = groupfemale - groupmale, 
                           levels=colnames(mm))
# Running contrast analysis
vfit <- contrasts.fit(fit, contrasts = contrasts)
# Looking at N of DEGs with adj. p <0.01 and log2FC>2
sumTable <- summary(decideTests(vfit, adjust.method = "BH", p.value = 0.05, lfc = 1))
sumTable
# Computing differential expression based on the empirical Bayes moderation of
# the standard errors towards a common value. Robust = should the estimation of
# the empirical Bayes prior parameters be robustified against outlier sample
# variances?
veBayesFit <- eBayes(vfit, robust=TRUE)
png(filename ="FIGURES/voomFits/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/decidua_voomFit_geneLevel_FinalModel_trans.png",width=6,height=6,units="in",res=1200)
plotSA(veBayesFit, main = "Final model: Mean-variance trend")
dev.off()
dev.off()

# Getting the DEGs. Log2FC of 1 is equivalent to linear fold change of 2.
# Getting summary statistics for all genes
allComparisons <- colnames(contrasts)
coef = 1

for (i in allComparisons){
  vTopTableAll <- topTable(veBayesFit, coef=coef, n=Inf, p.value=1, lfc=0)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/DEGs_deciduas_trans", i, "_fdr1_lfc0.txt", sep = "")
  write.table(vTopTableAll, path, sep = "\t")
  
  # Adj.p<0.05, log2FC>0
  vTopTableFDR.05 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.05, lfc=0)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/DEGs_deciduas_trans", i, "_fdr05_lfc0.txt", sep = "")
  write.table(vTopTableFDR.05, path, sep = "\t")
  
  # log2FC>1
  vTopTableFC1 <- topTable(veBayesFit, coef=coef, n=Inf, lfc=1)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/DEGs_deciduas_trans", i, "_lfc1.txt", sep = "")
  write.table(vTopTableFC1, path, sep = "\t")
  
  # log2FC>2
  vTopTableFC2 <- topTable(veBayesFit, coef=coef, n=Inf, lfc=2)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/DEGs_deciduas_trans", i, "_lfc2.txt", sep = "")
  write.table(vTopTableFC2, path, sep = "\t")
  
  # Adj.p<0.05, log2FC>1
  vTopTableFDR.05FC1 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.05, lfc=1)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/DEGs_deciduas_trans", i, "_fdr05_lfc1.txt", sep = "")
  write.table(vTopTableFDR.05FC1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>0
  vTopTableFDR.01 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=0)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/DEGs_deciduas_trans", i, "_fdr001_lfc0.txt", sep = "")
  write.table(vTopTableFDR.01, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>1
  vTopTableFDR.01FC1 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=1)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/DEGs_deciduas_trans", i, "_fdr001_lfc1.txt", sep = "")
  write.table(vTopTableFDR.01FC1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>2
  vTopTableFDR.01FC2 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=2)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/DEGs_deciduas_trans", i, "_fdr001_lfc2.txt", sep = "")
  write.table(vTopTableFDR.01FC2, path, sep = "\t")
  
  coef = coef+1
}

top.table <- topTable(veBayesFit, sort.by = "P", n = Inf)
head(top.table, 10)
length(which(top.table$adj.P.Val < 0.05))
summary(decideTests(veBayesFit))
dt <- decideTests(veBayesFit)

MvsF_All <- topTable(veBayesFit, n=Inf, coef=1)
P.Value <- c(MvsF_All$P.Value)
logFC <- c(MvsF_All$logFC)
adj.P.Val <- c(MvsF_All$adj.P.Val)
df <- data.frame(P.Value, logFC, adj.P.Val, MvsF_All$Geneid, MvsF_All$chr, MvsF_All$gene_name)

# Male bias
df.MB <- subset(df, adj.P.Val < 0.05 & logFC < -1 & MvsF_All$chr != "chrY") #define male-biased, black
df.MB <- cbind(df.MB, rep(3, nrow(df.MB)))
colnames(df.MB)[7] <- "color"

# not signigicant 
df.N <- subset(df, adj.P.Val >= 0.05) #define non-significant, grey
df.N <- cbind(df.N, rep(1, nrow(df.N)))
colnames(df.N)[7] <- "color"

df.S <- subset(df, (adj.P.Val < 0.05 & logFC > -1) | (adj.P.Val < 0.05 & logFC < 1)) #define non-significant, grey
df.S <- cbind(df.S, rep(2, nrow(df.S)))
colnames(df.S)[7] <- "color"

# Female bias
df.FB <- subset(df,  adj.P.Val < 0.05 & logFC > 1 & MvsF_All$chr != "chrX") #define female-biased, black
df.FB <- cbind(df.FB, rep(3, nrow(df.FB)))
colnames(df.FB)[7] <- "color"

# X-linked and significant 
df.X <- subset(df, (adj.P.Val < 0.05 & MvsF_All$chr == "chrX" & logFC > 1) | (adj.P.Val < 0.05 & MvsF_All$chr == "chrX" & logFC < -1)) #define X-chromosomal, red
df.X <- cbind(df.X, rep(4, nrow(df.X)))
colnames(df.X)[7] <- "color"

# Y-linked and signigicant 
df.Y <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr == "chrY" & logFC < -1) #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(5, nrow(df.Y)))
colnames(df.Y)[7] <- "color"

df.t <- rbind(df.N, df.S, df.MB, df.FB, df.X, df.Y)
df.t$color <- as.factor(df.t$color)
colnames(df.t)[4] <- "Geneid"
colnames(df.t)[5] <- "chr"
colnames(df.t)[6] <- "gene_name"
colnames(df.t)[7] <- "color"
dim(df.t)

dim(df)
dim(df.N)
dim(df.S)
dim(df.MB)
dim(df.Y)
dim(df.FB)
dim(df.X)

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 65)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") 
                      | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") 
                      | adj.P.Val<0.05 & logFC < -1 & (chr=="chrX") & (color=="4") 
                      |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))

#subset_data <- subset(df.t, (logFC < -1.5 & P.Value <0.01 & color=="1") | (logFC > 1.5 & P.Value <0.01 & color=="1"))


# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/deciduas_tranLevel_fdr05_lfc1_adjpsizeSet.png", width= 480, height= 580)
p + ggtitle("decidua differential expression transcript level")+
  # geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene_name), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  #  xlim(c(-10, 10)) + ylim(c(0, 40)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") 
                      | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") 
                      | adj.P.Val<0.05 & logFC < -1 & (chr=="chrX") & (color=="4") 
                      |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/deciduas_tranLevel_fdr05_lfc1.png", width= 480, height= 580)
p + ggtitle("decidua differential expression transcript level")+
  # geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene_name), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

MvsF_All <- topTable(veBayesFit, n=Inf, coef=1)
P.Value <- c(MvsF_All$P.Value)
logFC <- c(MvsF_All$logFC)
adj.P.Val <- c(MvsF_All$adj.P.Val)
df <- data.frame(P.Value, logFC, adj.P.Val, MvsF_All$Geneid, MvsF_All$chr, MvsF_All$gene_name)

# Male bias
df.MB <- subset(df, logFC < -1 & MvsF_All$chr != "chrY") #define male-biased, black
df.MB <- cbind(df.MB, rep(3, nrow(df.MB)))
colnames(df.MB)[7] <- "color"

# not signigicant 
df.N <- subset(df, adj.P.Val >= 0.05) #define non-significant, grey
df.N <- cbind(df.N, rep(1, nrow(df.N)))
colnames(df.N)[7] <- "color"

df.S <- subset(df, logFC > -1 | logFC < 1) #define non-significant, grey
df.S <- cbind(df.S, rep(2, nrow(df.S)))
colnames(df.S)[7] <- "color"

# Female bias
df.FB <- subset(df, logFC > 1 & MvsF_All$chr != "chrX") #define female-biased, black
df.FB <- cbind(df.FB, rep(3, nrow(df.FB)))
colnames(df.FB)[7] <- "color"

# X-linked and significant 
df.X <- subset(df, (MvsF_All$chr == "chrX" & logFC > 1) | (MvsF_All$chr == "chrX" & logFC < -1)) #define X-chromosomal, red
df.X <- cbind(df.X, rep(4, nrow(df.X)))
colnames(df.X)[7] <- "color"

# Y-linked and signigicant 
df.Y <- subset(df,  MvsF_All$chr == "chrY" & logFC < -1) #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(5, nrow(df.Y)))
colnames(df.Y)[7] <- "color"

df.t <- rbind(df.N, df.S, df.MB, df.FB, df.X, df.Y)
df.t$color <- as.factor(df.t$color)
colnames(df.t)[4] <- "Geneid"
colnames(df.t)[5] <- "chr"
colnames(df.t)[6] <- "gene_name"
colnames(df.t)[7] <- "color"
dim(df.t)

dim(df)
dim(df.N)
dim(df.S)
dim(df.MB)
dim(df.Y)
dim(df.FB)
dim(df.X)

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 65)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, logFC < -1 & (chr=="chrY") & (color=="5") 
                      |logFC > 1 & (chr=="chrX") & (color=="4") 
                      | logFC < -1 & (chr=="chrX") & (color=="4") 
                      | (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/deciduas_tranLevel_lfc1_sizeSet.png", width= 480, height= 580)
p + ggtitle("decidua differential expression transcript level")+
  # geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene_name), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  #  xlim(c(-10, 10)) + ylim(c(0, 40)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, logFC < -1 & (chr=="chrY") & (color=="5") 
                      | logFC > 1 & (chr=="chrX") & (color=="4") 
                      | logFC < -1 & (chr=="chrX") & (color=="4") 
                      | (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/deciduas_tranLevel_lfc1.png", width= 480, height= 580)
p + ggtitle("decidua differential expression transcript level")+
  # geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene_name), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

MvsF_All <- topTable(veBayesFit, n=Inf, coef=1)
P.Value <- c(MvsF_All$P.Value)
logFC <- c(MvsF_All$logFC)
adj.P.Val <- c(MvsF_All$adj.P.Val)
df <- data.frame(P.Value, logFC, adj.P.Val, MvsF_All$Geneid, MvsF_All$chr, MvsF_All$gene_name)

# Male bias
df.MB <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr != "chrY") #define male-biased, black
df.MB <- cbind(df.MB, rep(3, nrow(df.MB)))
colnames(df.MB)[7] <- "color"

# not signigicant 
df.N <- subset(df, adj.P.Val >= 0.05) #define non-significant, grey
df.N <- cbind(df.N, rep(1, nrow(df.N)))
colnames(df.N)[7] <- "color"

df.S <- subset(df, adj.P.Val < 0.05) #define non-significant, grey
df.S <- cbind(df.S, rep(2, nrow(df.S)))
colnames(df.S)[7] <- "color"

# Female bias
df.FB <- subset(df,  adj.P.Val < 0.05 & MvsF_All$chr != "chrX") #define female-biased, black
df.FB <- cbind(df.FB, rep(3, nrow(df.FB)))
colnames(df.FB)[7] <- "color"

# X-linked and significant 
df.X <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr == "chrX" ) #define X-chromosomal, red
df.X <- cbind(df.X, rep(4, nrow(df.X)))
colnames(df.X)[7] <- "color"

# Y-linked and signigicant 
df.Y <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr == "chrY") #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(5, nrow(df.Y)))
colnames(df.Y)[7] <- "color"

df.t <- rbind(df.N, df.S, df.MB, df.FB, df.X, df.Y)
df.t$color <- as.factor(df.t$color)
colnames(df.t)[4] <- "Geneid"
colnames(df.t)[5] <- "chr"
colnames(df.t)[6] <- "gene_name"
colnames(df.t)[7] <- "color"
dim(df.t)

dim(df)
dim(df.N)
dim(df.S)
dim(df.MB)
dim(df.Y)
dim(df.FB)
dim(df.X)

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 65)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, adj.P.Val<0.05 & (chr=="chrY") & (color=="5") 
                      | adj.P.Val<0.05 & (chr=="chrX") & (color=="4") 
                      | adj.P.Val<0.05 & (chr=="chrX") & (color=="4") 
                      |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/deciduas_tranLevel_fdr05_sizeSet.png", width= 480, height= 580)
p + ggtitle("decidua differential expression transcript level")+
  #  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene_name), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  #  xlim(c(-10, 10)) + ylim(c(0, 40)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, adj.P.Val<0.05  & (chr=="chrY") & (color=="5") 
                      | adj.P.Val<0.05  & (chr=="chrX") & (color=="4") 
                      | adj.P.Val<0.05  & (chr=="chrX") & (color=="4") 
                      |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA/deciduas_tranLevel_fdr05.png", width= 480, height= 580)
p + ggtitle("decidua differential expression transcript level")+
  #geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene_name), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------

#------------- VOOM 
#--------- gene level Voom 
group <-interaction(DGE_gene$samples$Offspring_sex) #, DGE$samples$batch, DGE$samples$Lane, DGE$samples$Reported_race) #, DGE$samples$batch, DGE$samples$Reported_race, DGE$samples$site, DGE$samples$Lane, DGE$samples)
#mm <- model.matrix(~0 + group + Maternal_age + Gravidity + Lane + PC1 + PC2)
mm <- model.matrix(~0 + group + Maternal_age + Gravidity + Gestational_age + Lane + PC1 + PC2 + MethodConception)
colnames(mm) <- gsub("v\\$targets\\$Offspring_sex", "", colnames(mm))
dim(mm)

y <- voom(DGE_gene, mm, plot = T) # plot T for plotting 

# Fitting linear models in limma
fit <- lmFit(y, mm)

png(filename ="FIGURES/voomFits/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/decidua_voomFit.png",width=6,height=6,units="in",res=1200)
v <- voomWithQualityWeights(DGE_gene, mm, plot=TRUE)
dev.off()
dev.off()

voomCounts <- v$E
voomCountsMatrix <- data.matrix(voomCounts, rownames.force = NA)
#write.table(voomCounts, "decidua_voom_cpm_genelevel.txt", sep="\t", row.names = FALSE, quote = FALSE)
genes <- v$genes
genesVoom <- cbind(genes, voomCounts)


#write.table(genes, "decidua_genes.txt", row.names = FALSE, quote = FALSE)
#v$genes
fit <- lmFit(v, mm)
png(filename ="FIGURES/voomFits/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/decidua_fitAmeanHist_fpkm1.png",width=6,height=6,units="in",res=1200)
hist(fit$Amean) 
dev.off()
dev.off()
plotSA(fit)
# keep <- fit$Amean > 0
# fit2 <- eBayes(fit[keep,], trend=TRUE) 
# hist(fit2$Amean)
# plotSA(fit2)

contrasts <- makeContrasts(FvsM = groupfemale - groupmale, 
                           levels=colnames(mm))
# Running contrast analysis
vfit <- contrasts.fit(fit, contrasts = contrasts)
# Looking at N of DEGs with adj. p <0.01 and log2FC>2
sumTable <- summary(decideTests(vfit, adjust.method = "BH", p.value = 0.05, lfc = 1))
sumTable
# Computing differential expression based on the empirical Bayes moderation of
# the standard errors towards a common value. Robust = should the estimation of
# the empirical Bayes prior parameters be robustified against outlier sample
# variances?
veBayesFit <- eBayes(vfit, robust=TRUE)
png(filename ="FIGURES/voomFits/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/decidua_voomFit_geneLevel_FinalModel.png",width=6,height=6,units="in",res=1200)
plotSA(veBayesFit, main = "Final model: Mean-variance trend")
dev.off()
dev.off()

# Getting the DEGs. Log2FC of 1 is equivalent to linear fold change of 2.
# Getting summary statistics for all genes
allComparisons <- colnames(contrasts)
coef = 1
#df <- data.frame(veBayesFit$sigma, veBayesFit$genes)
#df$sqrtSigma <- sqrt(df$veBayesFit.sigma)

for (i in allComparisons){
  vTopTableAll <- topTable(veBayesFit, coef=coef, n=Inf, p.value=1, lfc=0)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/DEGs_deciduas_genes", i, "_fdr1_lfc0.txt", sep = "")
  write.table(vTopTableAll, path, sep = "\t")
  
  # Adj.p<0.05, log2FC>0
  vTopTableFDR.05 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.05, lfc=0)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/DEGs_deciduas_genes", i, "_fdr05_lfc0.txt", sep = "")
  write.table(vTopTableFDR.05, path, sep = "\t")
  
  # log2FC>1
  vTopTableFC1 <- topTable(veBayesFit, coef=coef, n=Inf, lfc=1)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/DEGs_deciduas_genes", i, "_lfc1.txt", sep = "")
  write.table(vTopTableFC1, path, sep = "\t")
  
  # log2FC>2
  vTopTableFC2 <- topTable(veBayesFit, coef=coef, n=Inf, lfc=2)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/DEGs_deciduas_genes", i, "_lfc2.txt", sep = "")
  write.table(vTopTableFC2, path, sep = "\t")
  
  # Adj.p<0.05, log2FC>1
  vTopTableFDR.05FC1 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.05, lfc=1)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/DEGs_deciduas_genes", i, "_fdr05_lfc1.txt", sep = "")
  write.table(vTopTableFDR.05FC1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>0
  vTopTableFDR.01 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=0)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/DEGs_deciduas_genes", i, "_fdr001_lfc0.txt", sep = "")
  write.table(vTopTableFDR.01, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>1
  vTopTableFDR.01FC1 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=1)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/DEGs_deciduas_genes", i, "_fdr001_lfc1.txt", sep = "")
  write.table(vTopTableFDR.01FC1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>2
  vTopTableFDR.01FC2 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=2)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/DEGs_deciduas_genes", i, "_fdr001_lfc2.txt", sep = "")
  write.table(vTopTableFDR.01FC2, path, sep = "\t")
  
  coef = coef+1
}

top.table <- topTable(veBayesFit, sort.by = "P", n = Inf)
head(top.table, 10)
length(which(top.table$adj.P.Val < 0.05))
summary(decideTests(veBayesFit))
dt <- decideTests(veBayesFit)

MvsF_All <- topTable(veBayesFit, n=Inf, coef=1)
P.Value <- c(MvsF_All$P.Value)
logFC <- c(MvsF_All$logFC)
adj.P.Val <- c(MvsF_All$adj.P.Val)
df <- data.frame(P.Value, logFC, adj.P.Val, MvsF_All$Geneid, MvsF_All$chr)

# Male bias
df.MB <- subset(df, adj.P.Val < 0.05 & logFC < -1 & MvsF_All$chr != "chrY") #define male-biased, black
df.MB <- cbind(df.MB, rep(3, nrow(df.MB)))
colnames(df.MB)[6] <- "Color"

# not signigicant 
df.N <- subset(df, adj.P.Val >= 0.05) #define non-significant, grey
df.N <- cbind(df.N, rep(1, nrow(df.N)))
colnames(df.N)[6] <- "Color"

df.S <- subset(df, (adj.P.Val < 0.05 & logFC > -1) | (adj.P.Val < 0.05 & logFC < 1)) #define non-significant, grey
df.S <- cbind(df.S, rep(2, nrow(df.S)))
colnames(df.S)[6] <- "Color"

# Female bias
df.FB <- subset(df,  adj.P.Val < 0.05 & logFC > 1 & MvsF_All$chr != "chrX") #define female-biased, black
df.FB <- cbind(df.FB, rep(3, nrow(df.FB)))
colnames(df.FB)[6] <- "Color"

# X-linked and significant 
df.X <- subset(df, (adj.P.Val < 0.05 & MvsF_All$chr == "chrX" & logFC > 1) | (adj.P.Val < 0.05 & MvsF_All$chr == "chrX" & logFC < -1)) #define X-chromosomal, red
df.X <- cbind(df.X, rep(4, nrow(df.X)))
colnames(df.X)[6] <- "Color"

# Y-linked and signigicant 
df.Y <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr == "chrY" & logFC < -1) #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(5, nrow(df.Y)))
colnames(df.Y)[6] <- "Color"

df.t <- rbind(df.N, df.S, df.MB, df.FB, df.X, df.Y)
df.t$Color <- as.factor(df.t$Color)
colnames(df.t)[4] <- "Geneid"
colnames(df.t)[5] <- "chr"
colnames(df.t)[6] <- "color"
dim(df.t)

dim(df)
dim(df.N)
dim(df.S)
dim(df.MB)
dim(df.Y)
dim(df.FB)
dim(df.X)

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 65)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") 
                      | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") 
                      | adj.P.Val<0.05 & logFC < -1 & (chr=="chrX") & (color=="4") 
                      |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))

#subset_data <- subset(df.t, (logFC < -1.5 & P.Value <0.01 & color=="1") | (logFC > 1.5 & P.Value <0.01 & color=="1"))


# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/deciduas_geneLevel_fdr05_lfc1_adjpsizeSet.png", width= 480, height= 580)
p + ggtitle("decidua differential expression gene level")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  #  xlim(c(-10, 10)) + ylim(c(0, 40)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") 
                      | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") 
                      | adj.P.Val<0.05 & logFC < -1 & (chr=="chrX") & (color=="4") 
                      |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/deciduas_geneLevel_fdr05_lfc1.png", width= 480, height= 580)
p + ggtitle("decidua differential expression gene level")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

MvsF_All <- topTable(veBayesFit, n=Inf, coef=1)
P.Value <- c(MvsF_All$P.Value)
logFC <- c(MvsF_All$logFC)
adj.P.Val <- c(MvsF_All$adj.P.Val)
df <- data.frame(P.Value, logFC, adj.P.Val, MvsF_All$Geneid, MvsF_All$chr)

# Male bias
df.MB <- subset(df, logFC < -1 & MvsF_All$chr != "chrY") #define male-biased, black
df.MB <- cbind(df.MB, rep(3, nrow(df.MB)))
colnames(df.MB)[6] <- "Color"

# not signigicant 
df.N <- subset(df, adj.P.Val >= 0.05) #define non-significant, grey
df.N <- cbind(df.N, rep(1, nrow(df.N)))
colnames(df.N)[6] <- "Color"

df.S <- subset(df, logFC > -1 | logFC < 1) #define non-significant, grey
df.S <- cbind(df.S, rep(2, nrow(df.S)))
colnames(df.S)[6] <- "Color"

# Female bias
df.FB <- subset(df, logFC > 1 & MvsF_All$chr != "chrX") #define female-biased, black
df.FB <- cbind(df.FB, rep(3, nrow(df.FB)))
colnames(df.FB)[6] <- "Color"

# X-linked and significant 
df.X <- subset(df, (MvsF_All$chr == "chrX" & logFC > 1) | (MvsF_All$chr == "chrX" & logFC < -1)) #define X-chromosomal, red
df.X <- cbind(df.X, rep(4, nrow(df.X)))
colnames(df.X)[6] <- "Color"

# Y-linked and signigicant 
df.Y <- subset(df,  MvsF_All$chr == "chrY" & logFC < -1) #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(5, nrow(df.Y)))
colnames(df.Y)[6] <- "Color"

df.t <- rbind(df.N, df.S, df.MB, df.FB, df.X, df.Y)
df.t$Color <- as.factor(df.t$Color)
colnames(df.t)[4] <- "Geneid"
colnames(df.t)[5] <- "chr"
colnames(df.t)[6] <- "color"
dim(df.t)

dim(df)
dim(df.N)
dim(df.S)
dim(df.MB)
dim(df.Y)
dim(df.FB)
dim(df.X)

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 65)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, logFC < -1 & (chr=="chrY") & (color=="5") 
                      |logFC > 1 & (chr=="chrX") & (color=="4") 
                      | logFC < -1 & (chr=="chrX") & (color=="4") 
                      | (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/deciduas_geneLevel_lfc1_sizeSet.png", width= 480, height= 580)
p + ggtitle("decidua differential expression gene level")+
  #geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  #  xlim(c(-10, 10)) + ylim(c(0, 40)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, logFC < -1 & (chr=="chrY") & (color=="5") 
                      | logFC > 1 & (chr=="chrX") & (color=="4") 
                      | logFC < -1 & (chr=="chrX") & (color=="4") 
                      | (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/deciduas_geneLevel_lfc1.png", width= 480, height= 580)
p + ggtitle("decidua differential expression gene level")+
  #geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

MvsF_All <- topTable(veBayesFit, n=Inf, coef=1)
P.Value <- c(MvsF_All$P.Value)
logFC <- c(MvsF_All$logFC)
adj.P.Val <- c(MvsF_All$adj.P.Val)
df <- data.frame(P.Value, logFC, adj.P.Val, MvsF_All$Geneid, MvsF_All$chr)

# Male bias
df.MB <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr != "chrY") #define male-biased, black
df.MB <- cbind(df.MB, rep(3, nrow(df.MB)))
colnames(df.MB)[6] <- "Color"

# not signigicant 
df.N <- subset(df, adj.P.Val >= 0.05) #define non-significant, grey
df.N <- cbind(df.N, rep(1, nrow(df.N)))
colnames(df.N)[6] <- "Color"

df.S <- subset(df, adj.P.Val < 0.05) #define non-significant, grey
df.S <- cbind(df.S, rep(2, nrow(df.S)))
colnames(df.S)[6] <- "Color"

# Female bias
df.FB <- subset(df,  adj.P.Val < 0.05 & MvsF_All$chr != "chrX") #define female-biased, black
df.FB <- cbind(df.FB, rep(3, nrow(df.FB)))
colnames(df.FB)[6] <- "Color"

# X-linked and significant 
df.X <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr == "chrX" ) #define X-chromosomal, red
df.X <- cbind(df.X, rep(4, nrow(df.X)))
colnames(df.X)[6] <- "Color"

# Y-linked and signigicant 
df.Y <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr == "chrY") #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(5, nrow(df.Y)))
colnames(df.Y)[6] <- "Color"

df.t <- rbind(df.N, df.S, df.MB, df.FB, df.X, df.Y)
df.t$Color <- as.factor(df.t$Color)
colnames(df.t)[4] <- "Geneid"
colnames(df.t)[5] <- "chr"
colnames(df.t)[6] <- "color"
dim(df.t)

dim(df)
dim(df.N)
dim(df.S)
dim(df.MB)
dim(df.Y)
dim(df.FB)
dim(df.X)

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 65)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, adj.P.Val<0.05 & (chr=="chrY") & (color=="5") 
                      | adj.P.Val<0.05 & (chr=="chrX") & (color=="4") 
                      | adj.P.Val<0.05 & (chr=="chrX") & (color=="4") 
                      |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/deciduas_geneLevel_fdr05_sizeSet.png", width= 480, height= 580)
p + ggtitle("decidua differential expression gene level")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  #  xlim(c(-10, 10)) + ylim(c(0, 40)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, adj.P.Val<0.05  & (chr=="chrY") & (color=="5") 
                      | adj.P.Val<0.05  & (chr=="chrX") & (color=="4") 
                      | adj.P.Val<0.05  & (chr=="chrX") & (color=="4") 
                      |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/deciduas_geneLevel_fdr05.png", width= 480, height= 580)
p + ggtitle("decidua differential expression gene level")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

#------------- VOOM 
#--------- tran level Voom 
group <-interaction(DGE_tran$samples$Offspring_sex) #, DGE$samples$batch, DGE$samples$Lane, DGE$samples$Reported_race) #, DGE$samples$batch, DGE$samples$Reported_race, DGE$samples$site, DGE$samples$Lane, DGE$samples)
#mm <- model.matrix(~0 + group + Maternal_age + Gravidity + Lane + PC1 + PC2)
mm <- model.matrix(~0 + group + Maternal_age + Gravidity + Gestational_age + Lane + PC1 + PC2 + MethodConception)
colnames(mm) <- gsub("v\\$targets\\$Offspring_sex", "", colnames(mm))
dim(mm)

y <- voom(DGE_tran, mm, plot = T) # plot T for plotting 

# Fitting linear models in limma
fit <- lmFit(y, mm)

png(filename ="FIGURES/voomFits/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/decidua_voomFit_trans.png",width=6,height=6,units="in",res=1200)
v <- voomWithQualityWeights(DGE_tran, mm, plot=TRUE)
dev.off()
dev.off()

voomCounts <- v$E
voomCountsMatrix <- data.matrix(voomCounts, rownames.force = NA)
#write.table(voomCounts, "decidua_voom_cpm_genelevel.txt", sep="\t", row.names = FALSE, quote = FALSE)
genes <- v$genes
genesVoom <- cbind(genes, voomCounts)

#write.table(genes, "decidua_genes.txt", row.names = FALSE, quote = FALSE)
#v$genes
fit <- lmFit(v, mm)
png(filename ="FIGURES/voomFits/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/decidua_fitAmeanHist_fpkm1_trans.png",width=6,height=6,units="in",res=1200)
hist(fit$Amean) 
dev.off()
dev.off()
plotSA(fit)
# keep <- fit$Amean > 0
# fit2 <- eBayes(fit[keep,], trend=TRUE) 
# hist(fit2$Amean)
# plotSA(fit2)

contrasts <- makeContrasts(FvsM = groupfemale - groupmale, 
                           levels=colnames(mm))
# Running contrast analysis
vfit <- contrasts.fit(fit, contrasts = contrasts)
# Looking at N of DEGs with adj. p <0.01 and log2FC>2
sumTable <- summary(decideTests(vfit, adjust.method = "BH", p.value = 0.05, lfc = 1))
sumTable
# Computing differential expression based on the empirical Bayes moderation of
# the standard errors towards a common value. Robust = should the estimation of
# the empirical Bayes prior parameters be robustified against outlier sample
# variances?
veBayesFit <- eBayes(vfit, robust=TRUE)
png(filename ="FIGURES/voomFits/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/decidua_voomFit_geneLevel_FinalModel_trans.png",width=6,height=6,units="in",res=1200)
plotSA(veBayesFit, main = "Final model: Mean-variance trend")
dev.off()
dev.off()

# Getting the DEGs. Log2FC of 1 is equivalent to linear fold change of 2.
# Getting summary statistics for all genes
allComparisons <- colnames(contrasts)
coef = 1

for (i in allComparisons){
  vTopTableAll <- topTable(veBayesFit, coef=coef, n=Inf, p.value=1, lfc=0)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/DEGs_deciduas_trans", i, "_fdr1_lfc0.txt", sep = "")
  write.table(vTopTableAll, path, sep = "\t")
  
  # Adj.p<0.05, log2FC>0
  vTopTableFDR.05 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.05, lfc=0)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/DEGs_deciduas_trans", i, "_fdr05_lfc0.txt", sep = "")
  write.table(vTopTableFDR.05, path, sep = "\t")
  
  # log2FC>1
  vTopTableFC1 <- topTable(veBayesFit, coef=coef, n=Inf, lfc=1)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/DEGs_deciduas_trans", i, "_lfc1.txt", sep = "")
  write.table(vTopTableFC1, path, sep = "\t")
  
  # log2FC>2
  vTopTableFC2 <- topTable(veBayesFit, coef=coef, n=Inf, lfc=2)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/DEGs_deciduas_trans", i, "_lfc2.txt", sep = "")
  write.table(vTopTableFC2, path, sep = "\t")
  
  # Adj.p<0.05, log2FC>1
  vTopTableFDR.05FC1 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.05, lfc=1)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/DEGs_deciduas_trans", i, "_fdr05_lfc1.txt", sep = "")
  write.table(vTopTableFDR.05FC1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>0
  vTopTableFDR.01 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=0)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/DEGs_deciduas_trans", i, "_fdr001_lfc0.txt", sep = "")
  write.table(vTopTableFDR.01, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>1
  vTopTableFDR.01FC1 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=1)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/DEGs_deciduas_trans", i, "_fdr001_lfc1.txt", sep = "")
  write.table(vTopTableFDR.01FC1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>2
  vTopTableFDR.01FC2 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=2)
  path <- paste("./DEGs/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/DEGs_deciduas_trans", i, "_fdr001_lfc2.txt", sep = "")
  write.table(vTopTableFDR.01FC2, path, sep = "\t")
  
  coef = coef+1
}

top.table <- topTable(veBayesFit, sort.by = "P", n = Inf)
head(top.table, 10)
length(which(top.table$adj.P.Val < 0.05))
summary(decideTests(veBayesFit))
dt <- decideTests(veBayesFit)

MvsF_All <- topTable(veBayesFit, n=Inf, coef=1)
P.Value <- c(MvsF_All$P.Value)
logFC <- c(MvsF_All$logFC)
adj.P.Val <- c(MvsF_All$adj.P.Val)
df <- data.frame(P.Value, logFC, adj.P.Val, MvsF_All$Geneid, MvsF_All$chr, MvsF_All$gene_name)

# Male bias
df.MB <- subset(df, adj.P.Val < 0.05 & logFC < -1 & MvsF_All$chr != "chrY") #define male-biased, black
df.MB <- cbind(df.MB, rep(3, nrow(df.MB)))
colnames(df.MB)[7] <- "color"

# not signigicant 
df.N <- subset(df, adj.P.Val >= 0.05) #define non-significant, grey
df.N <- cbind(df.N, rep(1, nrow(df.N)))
colnames(df.N)[7] <- "color"

df.S <- subset(df, (adj.P.Val < 0.05 & logFC > -1) | (adj.P.Val < 0.05 & logFC < 1)) #define non-significant, grey
df.S <- cbind(df.S, rep(2, nrow(df.S)))
colnames(df.S)[7] <- "color"

# Female bias
df.FB <- subset(df,  adj.P.Val < 0.05 & logFC > 1 & MvsF_All$chr != "chrX") #define female-biased, black
df.FB <- cbind(df.FB, rep(3, nrow(df.FB)))
colnames(df.FB)[7] <- "color"

# X-linked and significant 
df.X <- subset(df, (adj.P.Val < 0.05 & MvsF_All$chr == "chrX" & logFC > 1) | (adj.P.Val < 0.05 & MvsF_All$chr == "chrX" & logFC < -1)) #define X-chromosomal, red
df.X <- cbind(df.X, rep(4, nrow(df.X)))
colnames(df.X)[7] <- "color"

# Y-linked and signigicant 
df.Y <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr == "chrY" & logFC < -1) #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(5, nrow(df.Y)))
colnames(df.Y)[7] <- "color"

df.t <- rbind(df.N, df.S, df.MB, df.FB, df.X, df.Y)
df.t$color <- as.factor(df.t$color)
colnames(df.t)[4] <- "Geneid"
colnames(df.t)[5] <- "chr"
colnames(df.t)[6] <- "gene_name"
colnames(df.t)[7] <- "color"
dim(df.t)

dim(df)
dim(df.N)
dim(df.S)
dim(df.MB)
dim(df.Y)
dim(df.FB)
dim(df.X)

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 65)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") 
                      | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") 
                      | adj.P.Val<0.05 & logFC < -1 & (chr=="chrX") & (color=="4") 
                      |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))

#subset_data <- subset(df.t, (logFC < -1.5 & P.Value <0.01 & color=="1") | (logFC > 1.5 & P.Value <0.01 & color=="1"))


# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/deciduas_tranLevel_fdr05_lfc1_adjpsizeSet.png", width= 480, height= 580)
p + ggtitle("decidua differential expression transcript level")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene_name), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  #  xlim(c(-10, 10)) + ylim(c(0, 40)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") 
                      | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") 
                      | adj.P.Val<0.05 & logFC < -1 & (chr=="chrX") & (color=="4") 
                      |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/deciduas_tranLevel_fdr05_lfc1.png", width= 480, height= 580)
p + ggtitle("decidua differential expression transcript level")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene_name), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

MvsF_All <- topTable(veBayesFit, n=Inf, coef=1)
P.Value <- c(MvsF_All$P.Value)
logFC <- c(MvsF_All$logFC)
adj.P.Val <- c(MvsF_All$adj.P.Val)
df <- data.frame(P.Value, logFC, adj.P.Val, MvsF_All$Geneid, MvsF_All$chr, MvsF_All$gene_name)

# Male bias
df.MB <- subset(df, logFC < -1 & MvsF_All$chr != "chrY") #define male-biased, black
df.MB <- cbind(df.MB, rep(3, nrow(df.MB)))
colnames(df.MB)[7] <- "color"

# not signigicant 
df.N <- subset(df, adj.P.Val >= 0.05) #define non-significant, grey
df.N <- cbind(df.N, rep(1, nrow(df.N)))
colnames(df.N)[7] <- "color"

df.S <- subset(df, logFC > -1 | logFC < 1) #define non-significant, grey
df.S <- cbind(df.S, rep(2, nrow(df.S)))
colnames(df.S)[7] <- "color"

# Female bias
df.FB <- subset(df, logFC > 1 & MvsF_All$chr != "chrX") #define female-biased, black
df.FB <- cbind(df.FB, rep(3, nrow(df.FB)))
colnames(df.FB)[7] <- "color"

# X-linked and significant 
df.X <- subset(df, (MvsF_All$chr == "chrX" & logFC > 1) | (MvsF_All$chr == "chrX" & logFC < -1)) #define X-chromosomal, red
df.X <- cbind(df.X, rep(4, nrow(df.X)))
colnames(df.X)[7] <- "color"

# Y-linked and signigicant 
df.Y <- subset(df,  MvsF_All$chr == "chrY" & logFC < -1) #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(5, nrow(df.Y)))
colnames(df.Y)[7] <- "color"

df.t <- rbind(df.N, df.S, df.MB, df.FB, df.X, df.Y)
df.t$color <- as.factor(df.t$color)
colnames(df.t)[4] <- "Geneid"
colnames(df.t)[5] <- "chr"
colnames(df.t)[6] <- "gene_name"
colnames(df.t)[7] <- "color"
dim(df.t)

dim(df)
dim(df.N)
dim(df.S)
dim(df.MB)
dim(df.Y)
dim(df.FB)
dim(df.X)

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 65)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, logFC < -1 & (chr=="chrY") & (color=="5") 
                      |logFC > 1 & (chr=="chrX") & (color=="4") 
                      | logFC < -1 & (chr=="chrX") & (color=="4") 
                      | (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/deciduas_tranLevel_lfc1_sizeSet.png", width= 480, height= 580)
p + ggtitle("decidua differential expression transcript level")+
  # geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene_name), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  #  xlim(c(-10, 10)) + ylim(c(0, 40)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, logFC < -1 & (chr=="chrY") & (color=="5") 
                      | logFC > 1 & (chr=="chrX") & (color=="4") 
                      | logFC < -1 & (chr=="chrX") & (color=="4") 
                      | (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/deciduas_tranLevel_lfc1.png", width= 480, height= 580)
p + ggtitle("decidua differential expression transcript level")+
  # geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene_name), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

MvsF_All <- topTable(veBayesFit, n=Inf, coef=1)
P.Value <- c(MvsF_All$P.Value)
logFC <- c(MvsF_All$logFC)
adj.P.Val <- c(MvsF_All$adj.P.Val)
df <- data.frame(P.Value, logFC, adj.P.Val, MvsF_All$Geneid, MvsF_All$chr, MvsF_All$gene_name)

# Male bias
df.MB <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr != "chrY") #define male-biased, black
df.MB <- cbind(df.MB, rep(3, nrow(df.MB)))
colnames(df.MB)[7] <- "color"

# not signigicant 
df.N <- subset(df, adj.P.Val >= 0.05) #define non-significant, grey
df.N <- cbind(df.N, rep(1, nrow(df.N)))
colnames(df.N)[7] <- "color"

df.S <- subset(df, adj.P.Val < 0.05) #define non-significant, grey
df.S <- cbind(df.S, rep(2, nrow(df.S)))
colnames(df.S)[7] <- "color"

# Female bias
df.FB <- subset(df,  adj.P.Val < 0.05 & MvsF_All$chr != "chrX") #define female-biased, black
df.FB <- cbind(df.FB, rep(3, nrow(df.FB)))
colnames(df.FB)[7] <- "color"

# X-linked and significant 
df.X <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr == "chrX" ) #define X-chromosomal, red
df.X <- cbind(df.X, rep(4, nrow(df.X)))
colnames(df.X)[7] <- "color"

# Y-linked and signigicant 
df.Y <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr == "chrY") #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(5, nrow(df.Y)))
colnames(df.Y)[7] <- "color"

df.t <- rbind(df.N, df.S, df.MB, df.FB, df.X, df.Y)
df.t$color <- as.factor(df.t$color)
colnames(df.t)[4] <- "Geneid"
colnames(df.t)[5] <- "chr"
colnames(df.t)[6] <- "gene_name"
colnames(df.t)[7] <- "color"
dim(df.t)

dim(df)
dim(df.N)
dim(df.S)
dim(df.MB)
dim(df.Y)
dim(df.FB)
dim(df.X)

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 65)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, adj.P.Val<0.05 & (chr=="chrY") & (color=="5") 
                      | adj.P.Val<0.05 & (chr=="chrX") & (color=="4") 
                      | adj.P.Val<0.05 & (chr=="chrX") & (color=="4") 
                      |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/deciduas_tranLevel_fdr05_sizeSet.png", width= 480, height= 580)
p + ggtitle("decidua differential expression transcript level")+
  #  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene_name), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha= 1, size = 2.5) +
  theme( legend.position = "none") +
  #  xlim(c(-10, 10)) + ylim(c(0, 40)) +
  scale_color_manual(values = c("azure3", "azure3", "#7570B3", "#66A61E", "#D95F02")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, adj.P.Val<0.05  & (chr=="chrY") & (color=="5") 
                      | adj.P.Val<0.05  & (chr=="chrX") & (color=="4") 
                      | adj.P.Val<0.05  & (chr=="chrX") & (color=="4") 
                      |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))

png(filename ="Figures/Volcano/deciduas/mat_Age_Gravity_Lane_PC1_PC2_GA_MethodConception/deciduas_tranLevel_fdr05.png", width= 480, height= 580)
p + ggtitle("decidua differential expression transcript level")+
  #geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene_name), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()









