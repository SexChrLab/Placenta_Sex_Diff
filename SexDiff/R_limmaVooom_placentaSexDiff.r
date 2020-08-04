#------------------------------------------------------------------------------------
# Single factor DE analysis using limma/voom
#
#       created by Kimberly Olney, July 19th 2019
#
#------------------------------------------------------------------------------------
#  
# To install packages use biocLite from bioconductor source
#source("http://bioconductor.org/biocLite.R")
#install.packages("BiocUpgrade")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Glimma")
#biocLite("mixOmics")
#biocLite("UncertainInterval")
#install.packages("edgeR")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("edgeR")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("variancePartition")
#biocLite("edgeR")
#biocLite("RColorBrewer")
#biocLite("dplyr")
#biocLite("ggplot2")
#biocLite("ggrepel")
#biocLite("devtools")
#biocLite("gplots") 
#biocLite("VennDiagram")
#library(limma)
#library(Glimma)
#setRepositories(cim)
#install.packages("devtools")
# install.packages("dplyr")
# source("https://www.bioconductor.org/biocLite.R")
# biocLite(c("alyssafrazee/RSkittleBrewer", "ballgown", "genefilter"))
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ballgown")

library(limma)
library(edgeR)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(mixOmics)
library(devtools)
library(dplyr)
library(RSkittleBrewer)
library(genefilter)
library(ballgown)
library(stringr)
library(ggrepel)
library(Glimma)
library(gplots)

#--------------
# set working directory 
#--------------
setwd("~/HISAT_FeatureCounts/batch1_and_batch2/")
#-----------
# Read in count, genes, and phenotype data
#-----------
counts <- read.delim("counts_pheno/batch1_and_batch2_transcriptCounts.tsv", header=TRUE, sep="\t")
colnames(counts) <- str_replace_all(colnames(counts), pattern="\\.","-") # replace . with - in sample names
#pheno <- read.csv("counts_pheno/batch1_and_batch2_FandM_pheno.csv", header=TRUE, sep=",")
pheno <- read.delim("counts_pheno/PLACENTA_PHENO_ALL.txt", header=TRUE, sep="\t")
#genes <- read.csv("genesID.csv", header=TRUE, sep = ",")
genes <- read.csv("counts_pheno/transcripts_gene_names.csv", header=TRUE, sep = "\t")
# batch 1 samples
batch1_sampleIDs <- c("OBG0044-1", "OBG0044-2", "OBG0053-1", "OBG0053-2", "OBG0068-1", "OBG0068-2", "OBG0111-1", "OBG0111-2", "OBG0112-1", "OBG0112-2", 
                     "OBG0115-1", "OBG0115-2", "OBG0116-1", "OBG0116-2", "OBG0117-1", "OBG0117-2", "OBG0118-1", "OBG0118-2", "OBG0120-1", "OBG0120-2", 
                     "OBG0122-1", "OBG0122-2", "OBG0123-1", "OBG0123-2", "OBG0126-1", "OBG0126-2", "OBG0130-1", "OBG0130-2", "OBG0132-1", "OBG0132-2", 
                     "OBG0133-1", "OBG0133-2", "OBG0156-1", "OBG0156-2", "OBG0158-1", "OBG0158-2", "OBG0166-1", "OBG0166-2", "OBG0170-1", "OBG0170-2", 
                     "OBG0174-1", "OBG0174-2", "OBG0175-1", "OBG0175-2", "OBG0178-1", "OBG0178-2", "YPOPS0006-1", "YPOPS0006-2")
# batch 2 samples
batch2_sampleIDs <- c("OBG0014-1", "OBG0014-2", "OBG0015-1", "OBG0015-2", "OBG0019-1", "OBG0019-2", "OBG0021-1", "OBG0021-2", "OBG0022-1", "OBG0022-2", "OBG0024-1", 
                      "OBG0024-2", "OBG0026-1", "OBG0026-2", "OBG0027-1", "OBG0027-2", "OBG0028-1", "OBG0028-2", "OBG0029-1", "OBG0029-2", "OBG0030-1", "OBG0030-2", 
                      "OBG0031-1", "OBG0031-2", "OBG0032-1", "OBG0032-2", "OBG0039-1", "OBG0039-2", "OBG0047-1", "OBG0047-2", "OBG0050-1", "OBG0050-2", "OBG0051-1", 
                      "OBG0051-2", "OBG0053B2-1", "OBG0053B2-2", "OBG0065-1", "OBG0065-2", "OBG0066-1", "OBG0066-2", "OBG0085-1", "OBG0085-2", "OBG0090-1", 
                      "OBG0090-2", "OBG0107-1", "OBG0107-2", "OBG0121-1", "OBG0121-2", "OBG0138-1", "OBG0138-2", "OBG0149-1", "OBG0149-2", "OBG0180-1", 
                      "OBG0180-2", "OBG0188-1", "OBG0188-2", "OBG0191-1", "OBG0191-2", "OBG0201-1", "OBG0201-2", "OBG0205-1", "OBG0205-2", "OBG0289-1", 
                      "OBG0289-2", "OBG0338-1", "OBG0338-2", "OBG0342-1", "OBG0342-2", "YPOPS0007M-1", "YPOPS0007M-2", "YPOPS0123M-1", "YPOPS0123M-2")
# samples to remove 
samplesRemove <- c("OBG0174-1", "OBG0174-2", "OBG0175-1", "OBG0175-2", "OBG0015-1", "OBG0015-2", "OBG0065-1", "OBG0065-2", "OBG0188-1", "OBG0188-2", 
                   "OBG0014-1", "OBG0014-2", "OBG0026-1","OBG0026-2", "YPOPS0007M-1","YPOPS0007M-2")

batch1_removals <- c("OBG0174-1", "OBG0174-2", "OBG0175-1", "OBG0175-2")
batch2_removals <- c("OBG0015-1", "OBG0015-2", "OBG0065-1", "OBG0065-2", "OBG0188-1", "OBG0188-2", "OBG0014-1", "OBG0014-2", "OBG0026-1","OBG0026-2", "YPOPS0007M-1","YPOPS0007M-2")

RNA_1_batch1 <- c("OBG0044-1", "OBG0053-1", "OBG0068-1","OBG0111-1", "OBG0112-1",
                     "OBG0115-1", "OBG0116-1", "OBG0117-1","OBG0118-1",  "OBG0120-1", 
                  "OBG0122-1", "OBG0123-1", "OBG0126-1", "OBG0130-1", "OBG0132-1",
                  "OBG0133-1", "OBG0156-1", "OBG0158-1", "OBG0166-1", "OBG0170-1",
                  "OBG0174-1", "OBG0175-1","OBG0178-1", "YPOPS0006-1" )
# batch 2 samples
RNA_1_batch2 <- c("OBG0014-1", "OBG0015-1", "OBG0019-1", "OBG0021-1", "OBG0022-1", "OBG0024-1", 
                  "OBG0026-1", "OBG0027-1", "OBG0028-1",  "OBG0029-1", "OBG0030-1", 
                  "OBG0031-1", "OBG0032-1", "OBG0039-1", "OBG0047-1", "OBG0050-1","OBG0051-1", 
                  "OBG0053B2-1","OBG0065-1","OBG0066-1", "OBG0085-1",  "OBG0090-1", 
                  "OBG0107-1", "OBG0121-1",  "OBG0138-1", "OBG0149-1", "OBG0180-1", 
                  "OBG0188-1",  "OBG0191-1",  "OBG0201-1", "OBG0205-1","OBG0289-1", 
                  "OBG0338-1", "OBG0342-1",  "YPOPS0007M-1", "YPOPS0123M-1")

RNA_2_batch1 <- c("OBG0044-2", "OBG0053-2", "OBG0068-2","OBG0111-2", "OBG0112-2",
                  "OBG0115-2", "OBG0116-2", "OBG0117-2","OBG0118-2",  "OBG0120-2", 
                  "OBG0122-2", "OBG0123-2", "OBG0126-2", "OBG0130-2", "OBG0132-2",
                  "OBG0133-2", "OBG0156-2", "OBG0158-2", "OBG0166-2", "OBG0170-2",
                  "OBG0174-2", "OBG0175-2","OBG0178-2", "YPOPS0006-2")
# batch 2 samples
RNA_2_batch2 <- c("OBG0014-2", "OBG0015-2", "OBG0019-2", "OBG0021-2", "OBG0022-2", "OBG0024-2", 
                  "OBG0026-2", "OBG0027-2", "OBG0028-2",  "OBG0029-2", "OBG0030-2", 
                  "OBG0031-2", "OBG0032-2", "OBG0039-2", "OBG0047-2", "OBG0050-2","OBG0051-2", 
                  "OBG0053B2-2","OBG0065-2","OBG0066-2", "OBG0085-2",  "OBG0090-2", 
                  "OBG0107-2", "OBG0121-2",  "OBG0138-2", "OBG0149-2", "OBG0180-2", 
                  "OBG0188-2",  "OBG0191-2",  "OBG0201-2", "OBG0205-2","OBG0289-2", 
                  "OBG0338-2", "OBG0342-2",  "YPOPS0007M-2", "YPOPS0123M-2")



samplesToRemove <- c(samplesRemove) # update depending on comparison being made
SAMPLE_LENGTH<-as.numeric(length(samplesToRemove)) # to call later 
half_sample_length <- SAMPLE_LENGTH/2 # half the sample length
# update counts data frame to exlude select samples 
removals<-(names(counts) %in% samplesToRemove[1:SAMPLE_LENGTH]) # for matching names create a value of TRUE or FALSES 
counts_ExRemovals <-counts[!removals] # create a new counts file that excludes (Ex) the removals IDs 

# update pheno data frame to exlude select samples 
pheno_ExRemovals <- pheno[! pheno$sample %in% samplesToRemove[1:SAMPLE_LENGTH],] # update 1:16 depending on size of samples to remove

# create DGE list object 
DGE <- DGEList(counts=counts_ExRemovals, genes = genes) # create DGEList object 
samplenames <- (pheno_ExRemovals$sample) # sample names are not unique because a sample may belong to multiple groups
colnames(DGE) <- samplenames

# create groups for the samples  

sex <- factor(pheno_ExRemovals$sex, levels=c("female", "male"))
batch <- factor(pheno_ExRemovals$batch, levels=c("1","2"))
race <- factor(pheno_ExRemovals$race, levels=c("Asian", "Black", "White", "Hispanic", "Other"))
site <- factor(pheno_ExRemovals$site, levels=c("RNA_1", "RNA_2"))
lane <- factor(pheno_ExRemovals$lane, levels=c("L001", "L002","L003", "L004","L005", "L006"))
rep <- factor(pheno_ExRemovals$REP, level=c("OBG0044", "OBG0053", "OBG0068", "OBG0111", "OBG0112", "OBG0115", "OBG0116", "OBG0117", "OBG0118", "OBG0120", "OBG0122", "OBG0123", "OBG0126", "OBG0130", "OBG0132", "OBG0133", "OBG0156", "OBG0158", "OBG0166", "OBG0170", "OBG0174", "OBG0175", "OBG0178", "YPOPS0006", "OBG0014", "OBG0015", "OBG0019", "OBG0021", "OBG0022", "OBG0024", "OBG0026", "OBG0027", "OBG0028", "OBG0029", "OBG0030", "OBG0031", "OBG0032", "OBG0039", "OBG0047", "OBG0050", "OBG0051", "OBG0053B2", "OBG0065", "OBG0066", "OBG0085", "OBG0090", "OBG0107", "OBG0121", "OBG0138", "OBG0149", "OBG0180", "OBG0188", "OBG0191", "OBG0201", "OBG0205", "OBG0289", "OBG0338", "OBG0342", "YPOPS0007M", "YPOPS0123M"))


# add groups to samples in dge list 
DGE$samples$sex <- sex
DGE$samples$batch <- batch
DGE$samples$race <- race
DGE$samples$site <- site
DGE$samples$lane <- lane
DGE$samples$rep <- rep
DGE$genes <- genes

# sum replciates 
#DGE <-sumTechReps(DGE, DGE$samples$rep) # comment out to not sum replicates 
sample_length_sumTech<-length(DGE$samples$group)
half_sample_length_sumTech <- sample_length_sumTech/2
DGE <- calcNormFactors(DGE) # Calculate normalization factors 
# NOTE calcNormFactors doesn???t normalize the data, it just calculates normalization factors for use downstream.
dim(DGE) 
DGEgeneNames<-unique(DGE$genes$gene_name)
#---- 
cpm <- cpm(DGE) #log2-transformed counts per gene per million mapped reads (cpm).
genesNames <- genes[1]
placenta_cpm <- cbind(genesNames, cpm)
write.table(placenta_cpm, "placenta_cpm.txt", sep ="\t", row.names = FALSE)


keep<- rowSums(cpm>1)>=half_sample_length_sumTech # in at least half the samples 
table(keep)
dge_CPM <- DGE[keep,, keep.lib.sizes=FALSE] # When you subset a DGEList and specify keep.lib.sizes=FALSE, the lib.size for each sample (cf. the y$samples data.frame) will be recalculated to be the sum of the counts left in the rows of the experiment for each sample.
dim(dge_CPM)
logdge_CPM <- cpm(dge_CPM, log=TRUE) # log transformting the data 
dim(logdge_CPM)

# Create a new variable ???group??? that combines pheno information
group <-interaction(dge_CPM$samples$sex) #, DGE$samples$batch, DGE$samples$lane, DGE$samples$race) #, DGE$samples$batch, DGE$samples$race, DGE$samples$site, DGE$samples$lane, DGE$samples)
batch <- dge_CPM$samples$batch

# Below specifies a model where each coefficient corresponds to a group mean
# Specify the model to be fitted. We do this before using voom since voom uses variances of the model residuals (observed - fitted)
# mm <- model.matrix(~0 + group) # To adjust for batch in the analysis, add batch to the end of the call to model matrix. Everything else about the code stays the same
mm <- model.matrix(~0 + group + batch) # for including batch in the model 
#mm <- model.matrix(~0 + group)
y <- voom(dge_CPM, mm, plot = F) # plot T for plotting 

# Fitting linear models in limma
fit <- lmFit(y, mm)

# Comparisons between groups (log fold-changes) are obtained as contrasts of these fitted linear models:
# Specify which groups to compare:
# contr <- makeContrasts(groupfemale - groupmale, levels = colnames(coef(fit)))
# F RNA_1 to M RNA_1
contr <- makeContrasts(groupfemale - groupmale, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr) # Estimate contrast for each gene or transcript
tmp <- eBayes(tmp) #Empirical Bayes smoothing of standard errors (shrinks standard errors that are much larger or smaller than those from other genes towards the average standard error)
plotSA(tmp, main="Final model: Mean-variance trend", ylim=c(0,3))

top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 10)
length(which(top.table$adj.P.Val < 0.05))
summary(decideTests(tmp))
dt <- decideTests(tmp)
FvsM_sig <-subset(top.table, adj.P.Val < 0.05)
# write.table(FvsM_sig, file="genelists/FvsM_sig_MergeReps_transcripts_batch1&2.txt",sep="\t", row.names=FALSE)
write.table(FvsM_sig, file="genelists/FvsM_sig_MergeReps_transcripts_batch1&2.txt",sep="\t", row.names=FALSE)

#glMDPlot(tmp, coef=1, status=dt, main=colnames(tmp)[1], counts=dge_CPM$counts, groups=group, launch=FALSE)

# Assigning color to highlight genes that have an absolute fold change > 2 
# and a p-value < Bonferroni cut-off
MvsF_All <- topTable(tmp, n=Inf, coef=1)
P.Value <- c(MvsF_All$P.Value)
logFC <- c(MvsF_All$logFC)
adj.P.Val <- c(MvsF_All$adj.P.Val)
df <- data.frame(P.Value, logFC, adj.P.Val, MvsF_All$gene_name, MvsF_All$chr)

# expression in the placenta regardless of sex. 
f_up_chrX<-subset(femaleUp, femaleUp$MvsF_All.chr == "chrX")
f_up_chrAuto <- subset(femaleUp, femaleUp$MvsF_All.chr != "chrX")
femaleUpgenes<-unique(femaleUp$MvsF_All.gene_name)

# expression in female placentas 
femaleUp <- (subset(df,df$logF > 0.5))
f_up_chrX<-subset(femaleUp, femaleUp$MvsF_All.chr == "chrX")
f_up_chrAuto <- subset(femaleUp, femaleUp$MvsF_All.chr != "chrX")
femaleUpgenes<-unique(femaleUp$MvsF_All.gene_name)

# expression in male placentas 
maleUp <- (subset(df,df$logF < -0.5))
m_up_chrY<-subset(maleUp, maleUp$MvsF_All.chr == "chrY")
m_up_chrX<-subset(maleUp, maleUp$MvsF_All.chr == "chrX")
m_up_chrAuto <- subset(maleUp, maleUp$MvsF_All.chr != "chrX" & maleUp$MvsF_All.chr != "chrY" )
maleUpgenes<-unique(maleUp$MvsF_All.gene_name)


df.MB <- subset(df, adj.P.Val < 0.05 & logFC > 1 & MvsF_All$chr != "chrY") #define male-biased, black
df.MB <- cbind(df.MB, rep(1, nrow(df.MB)))
colnames(df.MB)[6] <- "Color"

df.N <- subset(df, adj.P.Val >= 0.05) #define non-significant, grey
df.N<- cbind(df.N, rep(2, nrow(df.N)))
colnames(df.N)[6] <- "Color"

df.FB <- subset(df,  adj.P.Val < 0.05 & logFC < 1 & MvsF_All$chr != "chrX") #define female-biased, black
df.FB <- cbind(df.FB, rep(3, nrow(df.FB)))
colnames(df.FB)[6] <- "Color"

df.X <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr == "chrX") #define X-chromosomal, red
df.X <- cbind(df.X, rep(4, nrow(df.X)))
colnames(df.X)[6] <- "Color"

df.Y <- subset(df, adj.P.Val < 0.05 & MvsF_All$chr == "chrY") #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(5, nrow(df.Y)))
colnames(df.Y)[6] <- "Color"

df.t <- rbind(df.MB, df.N, df.FB, df.X, df.Y)
df.t$Color <- as.factor(df.t$Color)
colnames(df.t)[4] <- "gene_name"
colnames(df.t)[5] <- "chr"
colnames(df.t)[6] <- "color"

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=color ))+
  geom_point(alpha = 0.5, size = 1.75) +
  theme( legend.position = "none") +
 # xlim(c(-10, 10)) + ylim(c(0, 17)) +
  scale_color_manual(values = c("black", "azure3", "black", "red", "blue")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
p
dev.off()
dev.off()
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") | adj.P.Val<0.05 & logFC < -1 & (chr=="chrX") & (color=="4") |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3") & abs(logFC < -1))
# chr == "chr1" | chr=="chr2"| chr == "chr3" | chr=="chr4"| chr == "chr5" | chr=="chr6"| chr == "chr7" | chr=="chr8"| chr == "chr9" | chr=="chr10"| chr == "chr11" | chr=="chr12"| chr == "chr13" | chr=="chr14"| chr == "chr15" | chr=="chr16"| chr == "chr17" | chr=="chr18"| chr == "chr19" | chr=="chr20" | chr == "chr21" | chr=="chr22"
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
png(filename ="figures/MvsF_VOLCANOplotFC_MergeReps_transcripts_batch1&2.png", width= 480, height= 580)
p + ggtitle("Male and female placenta differential expression")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene_name))+geom_hline(yintercept = 2, colour="red", linetype="dashed"
) + geom_vline(xintercept = 1, colour="red", linetype="dashed"
) + geom_vline(xintercept = -1, colour="red", linetype="dashed")
dev.off()
dev.off()

# 
#plotMD(tmp, column=1, status=dt[,1], main=colnames(tmp)[1])#, xlim=c(-8,13))
# 
FvsM_topTable <- top.table$gene_name[1:20]
i <- which(y$genes$gene_name %in% FvsM_topTable)
mycol <- colorpanel(1000,"green","white","purple")

png(filename ="figures/MvsF_DE_Top20_heatmap_MergeReps_transcripts_batch1&2.png", width= 480, height= 580)
heatmap.2(logdge_CPM[i,], scale="row",
          labRow=y$genes$gene_name[i], labCol=y$targets$rep, 
          col=mycol, trace="none", density.info="none", 
          margin=c(6,6), lhei=c(2,10), dendrogram="column")
dev.off()
dev.off()

FvsM_topTable <- top.table$gene_name[1:100]
i <- which(y$genes$gene_name %in% FvsM_topTable)
mycol <- colorpanel(1000,"green","white","purple")

png(filename ="figures/MvsF_DE_Top100_heatmap_MergeReps_transcripts_batch1&2.png", width= 480, height= 580)
heatmap.2(logdge_CPM[i,], scale="row",
          labRow=y$genes$gene_name[i], labCol=y$targets$rep, 
          col=mycol, trace="none", density.info="none", 
          margin=c(6,6), lhei=c(2,10), dendrogram="column")
dev.off()
dev.off()

col.sex <- sex
levels(col.sex) <-  brewer.pal(nlevels(col.sex), "Set2")
col.sex <- as.character(col.sex)
png(filename ="figures/MDS_FandM_wholeTranscriptome_samples_cpm1_dim1&2_Top20_MergeReps_transcripts_batch1&2.png",width=6,height=6,units="in",res=1200)
plotMDS(DGE, labels = sex, top = 20, gene.selection = "common", col=col.sex, dim.plot = c(1,2))
title(main="MDS F and M whole transcriptome, PW Top 100")
dev.off()
dev.off()

png(filename ="figures/MDS_FandM_wholeTranscriptome_samples_cpm1_dim1&2_Top100_MergeReps_transcripts_batch1&2.png",width=6,height=6,units="in",res=1200)
plotMDS(logdge_CPM, labels = sex, top = 100, gene.selection = "common", col=col.sex, dim.plot = c(1,2))
title(main="MDS F and M whole transcriptome, PW Top 100")
dev.off()
dev.off()
