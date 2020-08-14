#   title: "placenta and decidua MDS plots"
# author: "Kimberly Olney"
# date: "02/04/2020"
# output:
#   pdf_document: default
# html_document:
#   df_print: paged

library(RColorBrewer)
library(matrixStats)
library(ggplot2)
library(edgeR)
library(DESeq)
library(limma)
library(doParallel)
library(variancePartition)
library(org.Hs.eg.db)
library(clusterProfiler)
library(GOSemSim)
library(biomaRt)
library(UpSetR)
library(VennDiagram)
library(ggrepel)
library(dplyr)
library(stringr)
library(forcats)

# set working directory 
setwd("~/Dropbox (ASU)/Placenta/BATCH2_PLACENTA_DECIDUA_ANALYSIS/HISAT_FeatureCounts/batch1_and_batch2/")

# Read in count, genes, and phenotype data
# For gene level and transcript level analysis
genes <- read.csv("counts_pheno/genesID.csv", header=TRUE, sep = ",") # gene ids
transcripts <- read.csv("counts_pheno/transcriptsID.csv", header=TRUE, sep = ",") # transcript ids

# placenta gene analysis information (donated by 'pt')
counts_pg <- read.delim("counts_pheno/pisarska_geneCounts.tsv", header=TRUE, sep="\t")
counts_pt <- read.delim("counts_pheno/pisarska_transcriptCounts.tsv", header=TRUE, sep="\t")
pheno_p <- read.csv("counts_pheno/pisarska_pheno.txt", header=TRUE, sep="\t")

no_removals <- c()

samplesToRemove <- c(no_removals) # update depending on comparison being made
SAMPLE_LENGTH <- as.numeric(length(samplesToRemove)) # to call later 
half_sample_length <- SAMPLE_LENGTH/2 # half the sample length

# update counts data frame to exlude select samples 
removals_pg <- (names(counts_pg) %in% samplesToRemove[1:SAMPLE_LENGTH]) # for matching names create a value of TRUE or FALSES 
counts_pg_ExRemovals <-counts_pg[!removals_pg] # create a new counts file that excludes (Ex) the removals IDs 

removals_pt <- (names(counts_pt) %in% samplesToRemove[1:SAMPLE_LENGTH]) # for matching names create a value of TRUE or FALSES 
counts_pt_ExRemovals <-counts_pt[!removals_pt] # create a new counts file that excludes (Ex) the removals IDs 

geneCounts <- cbind(genes, counts_pg_ExRemovals)
transcriptsCounts <- cbind(transcripts, counts_pt_ExRemovals)

# update pheno data frame to exlude select samples 
pheno_ExRemovals <- pheno_p[! pheno_p$sample %in% samplesToRemove[1:SAMPLE_LENGTH],] # update 1:16 depending on size of samples to remove

# create groups for the samples  
sex <- factor(pheno_ExRemovals$sex, levels=c("female", "male"))
sample <- factor(pheno_ExRemovals$sample)

# create DGE list object 
DGE_gene <- DGEList(counts=counts_pg_ExRemovals, genes = genes) # create DGEList object 
samplenames <- (pheno_ExRemovals$sample) # sample names are not unique because a sample may belong to multiple groups
colnames(DGE_gene) <- samplenames

DGE_tran <- DGEList(counts=counts_pt_ExRemovals, genes = transcripts) # create DGEList object 
samplenames <- (pheno_ExRemovals$sample) # sample names are not unique because a sample may belong to multiple groups
colnames(DGE_tran) <- samplenames

# add groups to samples in dge list 
DGE_gene$samples$sex <- sex
DGE_gene$samples$sample <- sample

# add groups to samples in dge list 
DGE_tran$samples$sex <- sex
DGE_tran$samples$sample <- sample

# data pre-processing
#---- 
DGE_gene <- calcNormFactors(DGE_gene, method="TMM")

fpkm <- rpkm(DGE_gene, gene.length=DGE_gene$genes$length)

female_mean_fpkm <- apply(as.data.frame(fpkm)
                          [(DGE_gene$samples$sex=="female")],
                          1, mean, na.rm=TRUE)
male_mean_fpkm <- apply(as.data.frame(fpkm)
                        [(DGE_gene$samples$sex=="male")],
                        1, mean, na.rm=TRUE)
# keep <- (female_mean_fpkm > 0.0 | male_mean_fpkm > 0.0)

keep <- (female_mean_fpkm > 1.0 | male_mean_fpkm > 1.0)
table(keep)
#------ 
DGE_gene <- DGE_gene[keep,,keep.lib.sizes=FALSE]
cpm_gene <- cpm(DGE_gene) #log2-transformed counts per gene per million mapped reads (cpm).
genesForCPM <- DGE_gene$genes
may <- cbind(genesForCPM, cpm_gene)
may$length <- NULL
may2 <- melt(may, id=c("Geneid", "chr"))
pheno_p <- subset(pheno_ExRemovals) #, site == "RNA_1")
pheno_p$variable <- pheno_p$sample
df_merged <- merge(may2, pheno_p, by.x="variable", all.x=TRUE)
write.table(df_merged, "genelists/pisarska/cpm_placenta_df_genes.txt", sep = "\t")
head(df_merged)

lcpm_gene <- cpm(cpm_gene, log=TRUE) 

keep<- rowSums(cpm_gene>0)>=0
table(keep)
dge_g_CPM <- DGE_gene[keep,, keep.lib.sizes=FALSE]
dim(dge_g_CPM)
logdge_g_CPM1 <- cpm(dge_g_CPM, log=TRUE) # log transformting the data 
dim(logdge_g_CPM1)

#---- 
DGE_tran <- calcNormFactors(DGE_tran, method="TMM")

fpkm <- rpkm(DGE_tran, gene.length=DGE_tran$genes$Length)

female_mean_fpkm <- apply(as.data.frame(fpkm)
                          [(DGE_tran$samples$sex=="female")],
                          1, mean, na.rm=TRUE)
male_mean_fpkm <- apply(as.data.frame(fpkm)
                        [(DGE_tran$samples$sex=="male")],
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


dim(logdge_g_CPM1)
#-------------------------------------------
col.sex <- sex
levels(col.sex) <-  brewer.pal(nlevels(col.sex), "Set2")
col.sex <- as.character(col.sex)
#-------------------------------------------

#--- Transcripts
png(filename ="FIGURES/MDS/pisarska/pisarska_transcriptCounts_common_top100transcripts_dim1&2.png",width=5,height=5,units="in",res=1200)
plotMDS(logdge_t_CPM1, labels=sex, gene.selection = "common", 
        top = 100, col=col.sex, dim.plot = c(1,2))
title(main="Pisarska
      Top 100 transcripts")
dev.off()
dev.off()


png(filename ="FIGURES/MDS/pisarska/pisarska_transcriptCounts_common_Alltranscripts_dim1&2.png",width=5,height=5,units="in",res=1200)
plotMDS(logdge_t_CPM1, labels=sex, gene.selection = "common", 
        col=col.sex, dim.plot = c(1,2))
title(main="Pisarska
      All transcripts")
dev.off()
dev.off()

#----------- genes 
png(filename ="FIGURES/MDS/pisarska/pisarska_geneCounts_common_Allgenes_dim1&2.png",width=5,height=5,units="in",res=1200)
plotMDS(logdge_g_CPM1, labels=sex, gene.selection = "common", 
        col=col.sex, dim.plot = c(1,2))
title(main="Pisarska
      All genes")
dev.off()
dev.off()

png(filename ="FIGURES/MDS/pisarska/pisarska_geneCounts_common_Top100genes_dim1&2.png",width=5,height=5,units="in",res=1200)
plotMDS(logdge_g_CPM1, labels=sex, gene.selection = "common",
        top = 100, col=col.sex, dim.plot = c(1,2))
title(main="Pisarska
      Top 100 genes")
dev.off()
dev.off()


#---- 2 by 2
png(filename ="FIGURES/MDS/pisarska/pisarska_dim1&2_PAR2x2.png",width=5,height=5,units="in",res=1200)
par(mfrow = c(2,2), mar = c(1,1,1,1) + 1.5)
plotMDS(logdge_g_CPM1, labels=sex, gene.selection = "common"
        , col=col.sex, dim.plot = c(1,2))
title(main="All genes")

plotMDS(logdge_g_CPM1, labels=sex, gene.selection = "common",
        top = 100, col=col.sex, dim.plot = c(1,2))
title(main="Top 100 genes")

plotMDS(logdge_t_CPM1, labels=sex, gene.selection = "common"
        , col=col.sex, dim.plot = c(1,2))
title(main="All transcripts")

plotMDS(logdge_t_CPM1, labels=sex, gene.selection = "common",
        top = 100, col=col.sex, dim.plot = c(1,2))
title(main="Top 100 transcripts")
dev.off()
dev.off()

#-------------------------------------------
# Voom
#-------------------------------------------
#--------- gene level Voom 
group <-interaction(dge_g_CPM$samples$sex) #, DGE$samples$batch, DGE$samples$lane, DGE$samples$race) #, DGE$samples$batch, DGE$samples$race, DGE$samples$site, DGE$samples$lane, DGE$samples)

mm <- model.matrix(~0 + group)
colnames(mm) <- gsub("v\\$targets\\$sex", "", colnames(mm))
head(mm)

y <- voom(dge_g_CPM, mm, plot = T) # plot T for plotting 

# Fitting linear models in limma
fit <- lmFit(y, mm)

png(filename ="FIGURES/voomFits/pisarska/placenta_voomFit_geneLevel.png",width=6,height=6,units="in",res=1200)
v <- voomWithQualityWeights(dge_g_CPM, mm, plot=TRUE)
dev.off()
dev.off()

voomCounts <- v$E
voomCountsMatrix <- data.matrix(voomCounts, rownames.force = NA)
genes <- v$genes
fit <- lmFit(v, mm)
contrasts <- makeContrasts(FvsM = groupfemale - groupmale, 
                           levels=colnames(mm))
# Running contrast analysis
vfit <- contrasts.fit(fit, contrasts = contrasts)
# Looking at N of DEGs with adj. p <0.01 and log2FC>2
sumTable <- summary(decideTests(vfit, adjust.method = "BH", p.value = 0.05, lfc = 1))

may <- cbind(genes, voomCounts)
may$length <- NULL
may2 <- melt(may, id=c("Geneid", "chr"))
df_merged <- merge(may2, pheno_ExRemovals, by.x="variable", by.y="sample", all.x=TRUE)
write.table(df_merged, "genelists/pisarska/vfit_placenta_df_genes.txt", sep = "\t")


# Computing differential expression based on the empirical Bayes moderation of
# the standard errors towards a common value. Robust = should the estimation of
# the empirical Bayes prior parameters be robustified against outlier sample
# variances?
veBayesFit <- eBayes(vfit, robust=TRUE)
png(filename ="FIGURES/voomFits/pisarska/placenta_voomFit_geneLevel_FinalModel.png",width=6,height=6,units="in",res=1200)
plotSA(veBayesFit, main = "Final model: Mean-variance trend")
dev.off()
dev.off()

# Getting the DEGs. Log2FC of 1 is equivalent to linear fold change of 2.
# Getting summary statistics for all genes
allComparisons <- colnames(contrasts)
coef = 1

# vTopTableAll <- topTable(veBayesFit, coef=1, n=Inf, p.value=1, lfc=0)
# path <- paste("./DEGs/pisarska/DEGs_placentas_genes_fdr1_lfc0_NOFPKMFILTER.txt", sep = "")
# write.table(vTopTableAll, path, sep = "\t")

for (i in allComparisons){
  vTopTableAll <- topTable(veBayesFit, coef=coef, n=Inf, p.value=1, lfc=0)
  path <- paste("./DEGs/pisarska/DEGs_placentas_genes", i, "_fdr1_lfc0.txt", sep = "")
  write.table(vTopTableAll, path, sep = "\t")
  
  # Adj.p<0.05, log2FC>0
  vTopTableFDR.05 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.05, lfc=0)
  path <- paste("./DEGs/pisarska/DEGs_placentas_genes", i, "_fdr05_lfc0.txt", sep = "")
  write.table(vTopTableFDR.05, path, sep = "\t")
  
  # log2FC>1
  vTopTableFC1 <- topTable(veBayesFit, coef=coef, n=Inf, lfc=1)
  path <- paste("./DEGs/pisarska/DEGs_placentas_genes", i, "_lfc1.txt", sep = "")
  write.table(vTopTableFC1, path, sep = "\t")
  
  # log2FC>2
  vTopTableFC2 <- topTable(veBayesFit, coef=coef, n=Inf, lfc=2)
  path <- paste("./DEGs/pisarska/DEGs_placentas_genes", i, "_lfc2.txt", sep = "")
  write.table(vTopTableFC2, path, sep = "\t")
  
  # Adj.p<0.05, log2FC>1
  vTopTableFDR.05FC1 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.05, lfc=1)
  path <- paste("./DEGs/pisarska/DEGs_placentas_genes", i, "_fdr05_lfc1.txt", sep = "")
  write.table(vTopTableFDR.05FC1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>0
  vTopTableFDR.01 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=0)
  path <- paste("./DEGs/pisarska/DEGs_placentas_genes", i, "_fdr001_lfc0.txt", sep = "")
  write.table(vTopTableFDR.01, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>1
  vTopTableFDR.01FC1 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=1)
  path <- paste("./DEGs/pisarska/DEGs_placentas_genes", i, "_fdr001_lfc1.txt", sep = "")
  write.table(vTopTableFDR.01FC1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>2
  vTopTableFDR.01FC2 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=2)
  path <- paste("./DEGs/pisarska/DEGs_placentas_genes", i, "_fdr001_lfc2.txt", sep = "")
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
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/pisarska/placenta_MvsF_geneLevel_fdr05_lfc1_adjpsizeSet.png", width= 480, height= 580)
p2 <- p + ggtitle("placenta differential expression gene level")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
p2
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
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/pisarska/placenta_MvsF_geneLevel_fdr05_lfc1.png", width= 480, height= 580)
p2 <- p + ggtitle("placenta differential expression gene level")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
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
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/pisarska/placenta_MvsF_geneLevel_lfc1_sizeSet.png", width= 480, height= 580)
p2 <- p + ggtitle("placenta differential expression gene level")+
 # geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
p2
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
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/pisarska/placenta_MvsF_geneLevel_lfc1.png", width= 480, height= 580)
p2 <- p + ggtitle("placenta differential expression gene level")+
 # geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
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
  geom_point(alpha= 1, size = 3.5) +
  theme( legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 65)) +
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
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/pisarska/placenta_MvsF_geneLevel_fdr05_sizeSet.png", width= 440, height= 540)
p2 <- p + ggtitle("")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), 
                  segment.alpha = 1, size = 5.5, fontface="italic") + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
p2
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
# subset_data <- subset(df.t, adj.P.Val<0.05  & (chr=="chrY") & (color=="5") 
#                       | adj.P.Val<0.05  & (chr=="chrX") & (color=="4") 
#                       | adj.P.Val<0.05  & (chr=="chrX") & (color=="4") 
#                       |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))
subset_data <- subset(df.t, logFC < -1 & (chr=="chrY") & (color=="5") 
                      | logFC > 1 & (chr=="chrX") & (color=="4") 
                      | logFC < -1 & (chr=="chrX") & (color=="4") 
                      | (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/pisarska/placenta_MvsF_geneLevel_fdr05.png", width= 480, height= 580)
p2 <- p + ggtitle("placenta differential expression gene level")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
p2
dev.off()
dev.off()


#--------- transcript level Voom 

#--------- transcript level Voom 
group <-interaction(dge_t_CPM$samples$sex) #, DGE$samples$batch, DGE$samples$lane, DGE$samples$race) #, DGE$samples$batch, DGE$samples$race, DGE$samples$site, DGE$samples$lane, DGE$samples)

mm <- model.matrix(~0 + group)
colnames(mm) <- gsub("v\\$targets\\$sex", "", colnames(mm))
head(mm)

y <- voom(dge_t_CPM, mm, plot = T) # plot T for plotting 

# Fitting linear models in limma
fit <- lmFit(y, mm)

png(filename ="FIGURES/voomFits/pisarska/placenta_voomFit_transLevel.png",width=6,height=6,units="in",res=1200)
v <- voomWithQualityWeights(dge_t_CPM, mm, plot=TRUE)
dev.off()
dev.off()

voomCounts <- v$E
voomCountsMatrix <- data.matrix(voomCounts, rownames.force = NA)
genes <- v$genes
fit <- lmFit(v, mm)
contrasts <- makeContrasts(FvsM = groupfemale - groupmale, 
                           levels=colnames(mm))
# Running contrast analysis
vfit <- contrasts.fit(fit, contrasts = contrasts)
# Looking at N of DEGs with adj. p <0.01 and log2FC>2
sumTable <- summary(decideTests(vfit, adjust.method = "BH", p.value = 0.05, lfc = 1))

# Computing differential expression based on the empirical Bayes moderation of
# the standard errors towards a common value. Robust = should the estimation of
# the empirical Bayes prior parameters be robustified against outlier sample
# variances?
veBayesFit <- eBayes(vfit, robust=TRUE)
png(filename ="FIGURES/voomFits/pisarska/placenta_voomFit_transLevel_FinalModel.png",width=6,height=6,units="in",res=1200)
plotSA(veBayesFit, main = "Final model: Mean-variance trend")
dev.off()
dev.off()

# Getting the DEGs. Log2FC of 1 is equivalent to linear fold change of 2.
# Getting summary statistics for all genes
allComparisons <- colnames(contrasts)
coef = 1

for (i in allComparisons){
  vTopTableAll <- topTable(veBayesFit, coef=coef, n=Inf, p.value=1, lfc=0)
  path <- paste("./DEGs/pisarska/DEGs_placentas_trans", i, "_fdr1_lfc0.txt", sep = "")
  write.table(vTopTableAll, path, sep = "\t")
  
  # Adj.p<0.05, log2FC>0
  vTopTableFDR.05 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.05, lfc=0)
  path <- paste("./DEGs/pisarska/DEGs_placentas_trans", i, "_fdr05_lfc0.txt", sep = "")
  write.table(vTopTableFDR.05, path, sep = "\t")
  
  # log2FC>1
  vTopTableFC1 <- topTable(veBayesFit, coef=coef, n=Inf, lfc=1)
  path <- paste("./DEGs/pisarska/DEGs_placentas_trans", i, "_lfc1.txt", sep = "")
  write.table(vTopTableFC1, path, sep = "\t")
  
  # log2FC>2
  vTopTableFC2 <- topTable(veBayesFit, coef=coef, n=Inf, lfc=2)
  path <- paste("./DEGs/pisarska/DEGs_placentas_trans", i, "_lfc2.txt", sep = "")
  write.table(vTopTableFC2, path, sep = "\t")
  
  # Adj.p<0.05, log2FC>1
  vTopTableFDR.05FC1 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.05, lfc=1)
  path <- paste("./DEGs/pisarska/DEGs_placentas_trans", i, "_fdr05_lfc1.txt", sep = "")
  write.table(vTopTableFDR.05FC1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>0
  vTopTableFDR.01 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=0)
  path <- paste("./DEGs/pisarska/DEGs_placentas_trans", i, "_fdr001_lfc0.txt", sep = "")
  write.table(vTopTableFDR.01, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>1
  vTopTableFDR.01FC1 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=1)
  path <- paste("./DEGs/pisarska/DEGs_placentas_trans", i, "_fdr001_lfc1.txt", sep = "")
  write.table(vTopTableFDR.01FC1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>2
  vTopTableFDR.01FC2 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=2)
  path <- paste("./DEGs/pisarska/DEGs_placentas_trans", i, "_fdr001_lfc2.txt", sep = "")
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
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/pisarska/placenta_MvsF_tranLevel_fdr05_lfc1_adjpsizeSet.png", width= 480, height= 580)
p + ggtitle("placenta differential expression gene level")+
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
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/pisarska/placenta_MvsF_tranLevel_fdr05_lfc1.png", width= 480, height= 580)
p + ggtitle("placenta differential expression gene level")+
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
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/pisarska/placenta_MvsF_tranLevel_lfc1_sizeSet.png", width= 480, height= 580)
p + ggtitle("placenta differential expression gene level")+
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
subset_data <- subset(df.t, logFC < -1 & (chr=="chrY") & (color=="5") 
                      | logFC > 1 & (chr=="chrX") & (color=="4") 
                      | logFC < -1 & (chr=="chrX") & (color=="4") 
                      | (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/pisarska/placenta_MvsF_tranLevel_lfc1.png", width= 480, height= 580)
p + ggtitle("placenta differential expression gene level")+
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
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/pisarska/placenta_MvsF_tranLevel_fdr05_sizeSet.png", width= 480, height= 580)
p + ggtitle("placenta differential expression gene level")+
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
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/pisarska/placenta_MvsF_tranLevel_fdr05.png", width= 480, height= 580)
p + ggtitle("placenta differential expression gene level")+
  #geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene_name), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()



