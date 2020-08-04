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

pheno_d <- read.csv("counts_pheno/decidua_pheno.csv", header=TRUE, sep=",")

# decidua samples to remove due to failed QC and/or outlier in MDS plot
decidua_removals <- c("YPOPS0007M-DEC", "OBG0021-DEC", "OBG0019-DEC")

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
offspring_sex <- factor(pheno_ExRemovals$offspring_sex, levels=c("female", "male"))
tissue <- factor(pheno_ExRemovals$tissue, levels=c("decidua", "placenta"))
batch <- factor(pheno_ExRemovals$batch, levels=c("1","2"))
race <- factor(pheno_ExRemovals$race, levels=c("Asian", "Black", "White", "Hispanic", "Other"))
site <- factor(pheno_ExRemovals$site, levels=c("RNA_1", "RNA_2"))
lane <- factor(pheno_ExRemovals$lane, levels=c("L001", "L002","L003", "L004","L005", "L006"))
PC1 <- factor(pheno_ExRemovals$PC1)
PC1 <- as.numeric(levels(PC1))[PC1]
PC2 <- factor(pheno_ExRemovals$PC2)
PC2 <- as.numeric(levels(PC2))[PC2]

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
DGE_gene$samples$offspring_sex <- offspring_sex
DGE_gene$samples$race <- race
DGE_gene$samples$site <- site
DGE_gene$samples$lane <- lane
DGE_gene$genes <- genes
DGE_gene$samples$PC1 <- PC1
DGE_gene$samples$PC2 <- PC2

# add groups to samples in dge list 
DGE_tran$samples$sex <- sex
DGE_tran$samples$sample <- sample
DGE_tran$samples$offspring_sex <- offspring_sex
DGE_tran$samples$race <- race
DGE_tran$samples$site <- site
DGE_tran$samples$lane <- lane
DGE_tran$genes <- genes
DGE_tran$samples$PC1 <- PC1
DGE_tran$samples$PC2 <- PC2

dim(DGE_gene) 
dim(DGE_tran) 

cpm_gene <- cpm(DGE_gene) #log2-transformed counts per gene per million mapped reads (cpm).
lcpm_gene <- cpm(cpm_gene, log=TRUE) 

cpm_tran <- cpm(DGE_tran) #log2-transformed counts per gene per million mapped reads (cpm).
lcpm_tran <- cpm(cpm_tran, log=TRUE) 

cpmBatch_tran <- removeBatchEffect(lcpm_tran, batch=DGE_tran$samples$PC2)


#-------------------------------------------
col.offspringSex <- offspring_sex
levels(col.offspringSex) <-  brewer.pal(nlevels(col.offspringSex), "Accent")
col.offspringSex <- as.character(col.offspringSex)

#-------------------------------------------
png(filename ="FIGURES/MDS/deciduas/deciduas_excludingFails_common_allgenes_dim1&2.png",width=5,height=5,units="in",res=1200)
plotMDS(lcpm_gene, labels=offspring_sex, gene.selection = "common", col=col.offspringSex, dim.plot = c(1,2))
title(main="All genes")
dev.off()
dev.off()

png(filename ="FIGURES/MDS/deciduas/deciduas_excludingFails_common_top100genes_dim1&2.png",width=5,height=5,units="in",res=1200)
plotMDS(lcpm_gene, labels=offspring_sex, top = 100, gene.selection = "common", col=col.offspringSex, dim.plot = c(1,2))
title(main="Top 100 genes")
dev.off()
dev.off()

png(filename ="FIGURES/MDS/deciduas/deciduas_excludingFails_common_allTranscripts_dim1&2.png",width=5,height=5,units="in",res=1200)
plotMDS(lcpm_tran, labels=offspring_sex, gene.selection = "common", col=col.offspringSex, dim.plot = c(1,2))
title(main="All transcripts")
dev.off()
dev.off()

png(filename ="FIGURES/MDS/deciduas/deciduas_excludingFails_common_top100transcripts_dim1&2.png",width=5,height=5,units="in",res=1200)
plotMDS(lcpm_tran, labels=offspring_sex, gene.selection = "common", top = 100, col=col.offspringSex, dim.plot = c(1,2))
title(main="Top 100 trasncripts")
dev.off()
dev.off()

png(filename ="FIGURES/MDS/deciduas/deciduas_excludingFails_dim1&2_PAR2x2.png",width=5,height=5,units="in",res=1200)
par(mfrow = c(2,2), mar = c(1,1,1,1) + 0.9)

plotMDS(lcpm_gene, labels=offspring_sex, gene.selection = "common", col=col.offspringSex, dim.plot = c(1,2))
title(main="All genes")

plotMDS(lcpm_gene, labels=offspring_sex, gene.selection = "common", top = 100, col=col.offspringSex, dim.plot = c(1,2))
title(main="Top 100 genes")

plotMDS(lcpm_tran, labels=offspring_sex, gene.selection = "common", col=col.offspringSex, dim.plot = c(1,2))
title(main="All transcripts")

plotMDS(lcpm_tran, labels=offspring_sex, gene.selection = "common", top = 100, col=col.offspringSex, dim.plot = c(1,2))
title(main="Top 100 transcripts")

dev.off()
dev.off()

#-------------------------------------------
# including failes 
#-------------------------------------------
# png(filename ="FIGURES/MDS/deciduas/deciduas__common_allgenes_dim1&2.png",width=5,height=5,units="in",res=1200)
# plotMDS(lcpm_gene, labels=offspring_sex, gene.selection = "common", col=col.offspringSex, dim.plot = c(1,2))
# title(main="All genes")
# dev.off()
# dev.off()
# 
# png(filename ="FIGURES/MDS/deciduas/deciduas__common_top100genes_dim1&2.png",width=5,height=5,units="in",res=1200)
# plotMDS(lcpm_gene, labels=offspring_sex, top = 100, gene.selection = "common", col=col.offspringSex, dim.plot = c(1,2))
# title(main="Top 100 genes")
# dev.off()
# dev.off()
# 
# png(filename ="FIGURES/MDS/deciduas/deciduas__common_allTranscripts_dim1&2.png",width=5,height=5,units="in",res=1200)
# plotMDS(lcpm_tran, labels=offspring_sex, gene.selection = "common", col=col.offspringSex, dim.plot = c(1,2))
# title(main="All transcripts")
# dev.off()
# dev.off()
# 
# png(filename ="FIGURES/MDS/deciduas/deciduas__common_top100transcripts_dim1&2.png",width=5,height=5,units="in",res=1200)
# plotMDS(lcpm_tran, labels=offspring_sex, gene.selection = "common", top = 100, col=col.offspringSex, dim.plot = c(1,2))
# title(main="Top 100 trasncripts")
# dev.off()
# dev.off()
# 
# png(filename ="FIGURES/MDS/deciduas/deciduas__dim1&2_PAR2x2.png",width=5,height=5,units="in",res=1200)
# par(mfrow = c(2,2), mar = c(1,1,1,1) + 0.9)
# 
# plotMDS(lcpm_gene, labels=offspring_sex, gene.selection = "common", col=col.offspringSex, dim.plot = c(1,2))
# title(main="All genes")
# 
# plotMDS(lcpm_gene, labels=offspring_sex, gene.selection = "common", top = 100, col=col.offspringSex, dim.plot = c(1,2))
# title(main="Top 100 genes")
# 
# plotMDS(lcpm_tran, labels=offspring_sex, gene.selection = "common", col=col.offspringSex, dim.plot = c(1,2))
# title(main="All transcripts")
# 
# plotMDS(lcpm_tran, labels=offspring_sex, gene.selection = "common", top = 100, col=col.offspringSex, dim.plot = c(1,2))
# title(main="Top 100 transcripts")
# 
# dev.off()
# dev.off()

