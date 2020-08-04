#   title: "placenta and decidua MDS plots"
# author: "Kimberly Olney"
# date: "02/04/2020"
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

# placenta gene analysis information (donated by 'pg')
# placenta gene analysis information (donated by 'pt')
counts_pg <- read.delim("counts_pheno/batch1_and_batch2_geneCounts.tsv", header=TRUE, sep="\t")
colnames(counts_pg) <- str_replace_all(colnames(counts_pg), pattern="\\.","-") # replace . with - in sample names
counts_pt <- read.delim("counts_pheno/batch1_and_batch2_transcriptCounts.tsv", header=TRUE, sep="\t")
colnames(counts_pt) <- str_replace_all(colnames(counts_pt), pattern="\\.","-") # replace . with - in sample names
pheno_p <- read.csv("counts_pheno/batch1_and_batch2_FandM_pheno.txt", header=TRUE, sep="\t")

# decidua gene analysis information (donated by 'dg')
counts_dg <- read.delim("counts_pheno/decidua_geneCounts.tsv", header=TRUE, sep="\t")
colnames(counts_dg) <- str_replace_all(colnames(counts_dg), pattern="\\.","-") # replace . with - in sample names
counts_dt <- read.delim("counts_pheno/decidua_transcriptCounts.tsv", header=TRUE, sep="\t")
colnames(counts_dt) <- str_replace_all(colnames(counts_dt), pattern="\\.","-") # replace . with - in sample names
pheno_d <- read.csv("counts_pheno/decidua_pheno.csv", header=TRUE, sep=",")

# samples in batch and outlier samples

# batch1 samples 
placenta_batch1_sampleIDs <- c("OBG0044-1", "OBG0044-2", "OBG0053-1", "OBG0053-2", "OBG0068-1", 
                       "OBG0068-2", "OBG0111-1", "OBG0111-2", "OBG0112-1", "OBG0112-2",
                       "OBG0115-1", "OBG0115-2", "OBG0116-1", "OBG0116-2", "OBG0117-1", 
                       "OBG0117-2", "OBG0118-1", "OBG0118-2", "OBG0120-1", "OBG0120-2", 
                       "OBG0122-1", "OBG0122-2", "OBG0123-1", "OBG0123-2", "OBG0126-1", 
                       "OBG0126-2", "OBG0130-1", "OBG0130-2", "OBG0132-1", "OBG0132-2", 
                       "OBG0133-1", "OBG0133-2", "OBG0156-1", "OBG0156-2", "OBG0158-1", 
                       "OBG0158-2", "OBG0166-1", "OBG0166-2", "OBG0170-1", "OBG0170-2", 
                       "OBG0174-1", "OBG0174-2", "OBG0175-1", "OBG0175-2", "OBG0178-1", 
                       "OBG0178-2", "YPOPS0006-1", "YPOPS0006-2")

# batch 2 samples
placenta_batch2_sampleIDs <- c("OBG0014-1", "OBG0014-2", "OBG0015-1", "OBG0015-2", "OBG0019-1", 
                      "OBG0019-2", "OBG0021-1", "OBG0021-2", "OBG0022-1", "OBG0022-2", 
                      "OBG0024-1", "OBG0024-2", "OBG0026-1", "OBG0026-2", "OBG0027-1", 
                      "OBG0027-2", "OBG0028-1", "OBG0028-2", "OBG0029-1", "OBG0029-2", 
                      "OBG0030-1", "OBG0030-2", "OBG0031-1", "OBG0031-2", "OBG0032-1", 
                      "OBG0032-2", "OBG0039-1", "OBG0039-2", "OBG0047-1", "OBG0047-2", 
                      "OBG0050-1", "OBG0050-2", "OBG0051-1", "OBG0051-2", "OBG0053B2-1", 
                      "OBG0053B2-2", "OBG0065-1", "OBG0065-2", "OBG0066-1", "OBG0066-2", 
                      "OBG0085-1", "OBG0085-2", "OBG0090-1","OBG0090-2", "OBG0107-1", 
                      "OBG0107-2", "OBG0121-1", "OBG0121-2", "OBG0138-1", "OBG0138-2", 
                      "OBG0149-1", "OBG0149-2", "OBG0180-1", "OBG0180-2", "OBG0188-1", 
                      "OBG0188-2", "OBG0191-1", "OBG0191-2", "OBG0201-1", "OBG0201-2", 
                      "OBG0205-1", "OBG0205-2", "OBG0289-1", "OBG0289-2", "OBG0338-1", 
                      "OBG0338-2", "OBG0342-1", "OBG0342-2", "YPOPS0007M-1", "YPOPS0007M-2", 
                      "YPOPS0123M-1", "YPOPS0123M-2")

# samples to remove due to failed QC and/or outlier in MDS plot
placenta_batch1_removals <- c("OBG0174-1", "OBG0174-2", "OBG0175-1", "OBG0175-2")
placenta_batch2_removals <- c("OBG0015-1", "OBG0015-2", "OBG0065-1", "OBG0065-2", 
                              "OBG0188-1", "OBG0188-2", "OBG0014-1", "OBG0014-2", 
                              "OBG0026-1","OBG0026-2", "YPOPS0007M-1","YPOPS0007M-2")

placenta_removals <- c(placenta_batch1_removals, placenta_batch2_removals)

# decidua samples to remove due to failed QC and/or outlier in MDS plot
decidua_removals <- c("YPOPS0007M-DEC", "OBG0021-DEC", "OBG0019-DEC")
# corresponding placenta samples to remove 
decidua_placenta_removals <- c("OBG0021-1", "OBG0021-2", "OBG0019-1", "OBG0019-2")


no_removals <- c()
all_removals <- c(placenta_batch1_sampleIDs, placenta_batch2_removals, 
                  decidua_removals, decidua_placenta_removals)

samplesToRemove <- c(no_removals) # update depending on comparison being made
SAMPLE_LENGTH <- as.numeric(length(samplesToRemove)) # to call later 
half_sample_length <- SAMPLE_LENGTH/2 # half the sample length

# placneta and decidua gene counts 'pdg'
# placneta and decidua transcript counts 'pdt'
# placenta and decidua pheno data 'pd'
counts_pdg <- cbind(counts_pg, counts_dg)
counts_pdt <- cbind(counts_pt, counts_dt)
pheno_pd <- rbind.fill(pheno_p, pheno_d)

# update counts data frame to exlude select samples 
removals_pdg <- (names(counts_pdg) %in% samplesToRemove[1:SAMPLE_LENGTH]) # for matching names create a value of TRUE or FALSES 
counts_pdg_ExRemovals <-counts_pdg[!removals_pdg] # create a new counts file that excludes (Ex) the removals IDs 

removals_pdt <- (names(counts_pdt) %in% samplesToRemove[1:SAMPLE_LENGTH]) # for matching names create a value of TRUE or FALSES 
counts_pdt_ExRemovals <-counts_pdt[!removals_pdt] # create a new counts file that excludes (Ex) the removals IDs 

geneCounts <- cbind(genes, counts_pdg_ExRemovals)
transcriptsCounts <- cbind(transcripts, counts_pdt_ExRemovals)

# update pheno data frame to exlude select samples 
pheno_ExRemovals <- pheno_pd[! pheno_pd$sample %in% samplesToRemove[1:SAMPLE_LENGTH],] # update 1:16 depending on size of samples to remove

# create groups for the samples  
sex <- factor(pheno_ExRemovals$sex, levels=c("female", "male"))
sample <- factor(pheno_ExRemovals$sample)
offspring_sex <- factor(pheno_ExRemovals$offspring_sex, levels=c("female", "male"))
tissue <- factor(pheno_ExRemovals$tissue, levels=c("decidua", "placenta"))
batch <- factor(pheno_ExRemovals$batch, levels=c("1","2"))
race <- factor(pheno_ExRemovals$race, levels=c("Asian", "Black", "White", "Hispanic", "Other"))
site <- factor(pheno_ExRemovals$site, levels=c("RNA_1", "RNA_2"))
lane <- factor(pheno_ExRemovals$lane, levels=c("L001", "L002","L003", "L004","L005", "L006"))
rep <- factor(pheno_ExRemovals$REP, level=c("OBG0044", "OBG0053", "OBG0068", "OBG0111", "OBG0112", "OBG0115", "OBG0116", "OBG0117", "OBG0118", "OBG0120", "OBG0122", "OBG0123", "OBG0126", "OBG0130", "OBG0132", "OBG0133", "OBG0156", "OBG0158", "OBG0166", "OBG0170", "OBG0174", "OBG0175", "OBG0178", "YPOPS0006", "OBG0014", "OBG0015", "OBG0019", "OBG0021", "OBG0022", "OBG0024", "OBG0026", "OBG0027", "OBG0028", "OBG0029", "OBG0030", "OBG0031", "OBG0032", "OBG0039", "OBG0047", "OBG0050", "OBG0051", "OBG0053B2", "OBG0065", "OBG0066", "OBG0085", "OBG0090", "OBG0107", "OBG0121", "OBG0138", "OBG0149", "OBG0180", "OBG0188", "OBG0191", "OBG0201", "OBG0205", "OBG0289", "OBG0338", "OBG0342", "YPOPS0007M", "YPOPS0123M"))
PC1 <- factor(pheno_ExRemovals$PC1)
PC1 <- as.numeric(levels(PC1))[PC1]
PC2 <- factor(pheno_ExRemovals$PC2)
PC2 <- as.numeric(levels(PC2))[PC2]

# create DGE list object 
DGE_gene <- DGEList(counts=counts_pdg_ExRemovals, genes = genes) # create DGEList object 
samplenames <- (pheno_ExRemovals$sample) # sample names are not unique because a sample may belong to multiple groups
colnames(DGE_gene) <- samplenames

DGE_tran <- DGEList(counts=counts_pdt_ExRemovals, genes = transcripts) # create DGEList object 
samplenames <- (pheno_ExRemovals$sample) # sample names are not unique because a sample may belong to multiple groups
colnames(DGE_tran) <- samplenames

# add groups to samples in dge list 
DGE_gene$samples$sex <- sex
DGE_gene$samples$sample <- sample
DGE_gene$samples$offspring_sex <- offspring_sex
DGE_gene$samples$batch <- batch
DGE_gene$samples$race <- race
DGE_gene$samples$site <- site
DGE_gene$samples$lane <- lane
DGE_gene$samples$rep <- rep
DGE_gene$genes <- genes
DGE_gene$samples$PC1 <- PC1
DGE_gene$samples$PC2 <- PC2

# add groups to samples in dge list 
DGE_tran$samples$sex <- sex
DGE_tran$samples$sample <- sample
DGE_tran$samples$offspring_sex <- offspring_sex
DGE_tran$samples$batch <- batch
DGE_tran$samples$race <- race
DGE_tran$samples$site <- site
DGE_tran$samples$lane <- lane
DGE_tran$samples$rep <- rep
DGE_tran$genes <- genes
DGE_tran$samples$PC1 <- PC1
DGE_tran$samples$PC2 <- PC2

# sum replciates 
#DGE <-sumTechReps(DGE, DGE$samples$rep) # comment out to not sum replicates 
#sample_length_sumTech<-length(DGE$samples$group)
#half_sample_length_sumTech <- sample_length_sumTech/2
#DGE <- calcNormFactors(DGE) # Calculate normalization factors 
# NOTE calcNormFactors doesn’t normalize the data, it just calculates normalization factors for use downstream.

dim(DGE_gene) 
dim(DGE_tran) 

cpm_gene <- cpm(DGE_gene) #log2-transformed counts per gene per million mapped reads (cpm).
lcpm_gene <- cpm(cpm_gene, log=TRUE) 

cpm_tran <- cpm(DGE_tran) #log2-transformed counts per gene per million mapped reads (cpm).
lcpm_tran <- cpm(cpm_tran, log=TRUE) 

# fpkm <- rpkm(DGE, gene.length=DGE$genes$Length)
# 
# female_mean_fpkm <- apply(as.data.frame(fpkm)
#                           [(DGE$samples$offspring_sex=="female")],
#                           1, mean, na.rm=TRUE)
# male_mean_fpkm <- apply(as.data.frame(fpkm)
#                         [(DGE$samples$offspring_sex=="male")],
#                         1, mean, na.rm=TRUE)
# 
# 
# keep <- (female_mean_fpkm > 0.0 | male_mean_fpkm > 0.0)
# DGE <- DGE[keep,,keep.lib.sizes=FALSE]
# DGE <- calcNormFactors(DGE, method="TMM")
# keep <- rowSums(DGE$counts > 6) >= 2
# DGE <- DGE[keep,,keep.lib.size=FALSE]
# DGE <- calcNormFactors(DGE, method="TMM")
# 
# dim(dge$genes) # N of genes retained after filtering

#-------------------------------------------
col.tissue <- tissue
levels(col.tissue) <-  brewer.pal(nlevels(col.tissue), "Set1")
col.tissue <- as.character(col.tissue)

#-------------------------------------------

png(filename ="FIGURES/MDS/placentas_deciduas_common_allgenes_dim1&2.png",width=5,height=5,units="in",res=1200)
plotMDS(lcpm_gene, labels=tissue, gene.selection = "common", col=col.tissue, dim.plot = c(1,2))
title(main="placentas and deciduas
      common gene selection
      all genes")
dev.off()
dev.off()

png(filename ="FIGURES/MDS/placentas_deciduas_common_top100genes_dim1&2.png",width=5,height=5,units="in",res=1200)
plotMDS(lcpm_gene, labels=tissue, gene.selection = "common", top = 100, col=col.tissue, dim.plot = c(1,2))
title(main="placentas and deciduas
      common gene selection
      top 100 gene")
dev.off()
dev.off()

png(filename ="FIGURES/MDS/placentas_deciduas_common_allTranscripts_dim1&2.png",width=5,height=5,units="in",res=1200)
plotMDS(lcpm_tran, labels=tissue, gene.selection = "common", col=col.tissue, dim.plot = c(1,2))
title(main="placentas and deciduas
      common gene selection
      all transcripts")
dev.off()
dev.off()

png(filename ="FIGURES/MDS/placentas_deciduas_common_top100transcripts_dim1&2.png",width=5,height=5,units="in",res=1200)
plotMDS(lcpm_tran, labels=tissue, gene.selection = "common", top = 100, col=col.tissue, dim.plot = c(1,2))
title(main="placentas and deciduas
      common gene selection
      top 100 transcripts")
dev.off()
dev.off()


# ----- excluding bad samples
# png(filename ="FIGURES/MDS/placentas_deciduas_excludeRemovals_common_allgenes_dim1&2.png",width=5,height=5,units="in",res=1200)
# plotMDS(lcpm_gene, labels=tissue, gene.selection = "common", col=col.tissue, dim.plot = c(1,2))
# title(main="placentas and deciduas
#       common gene selection
#       all genes")
# dev.off()
# dev.off()
# 
# png(filename ="FIGURES/MDS/placentas_deciduas_excludeRemovals_common_top100genes_dim1&2.png",width=5,height=5,units="in",res=1200)
# plotMDS(lcpm_gene, labels=tissue, gene.selection = "common", top = 100, col=col.tissue, dim.plot = c(1,2))
# title(main="placentas and deciduas
#       common gene selection
#       top 100 gene")
# dev.off()
# dev.off()
# 
# png(filename ="FIGURES/MDS/placentas_deciduas_excludeRemovals_common_allTranscripts_dim1&2.png",width=5,height=5,units="in",res=1200)
# plotMDS(lcpm_tran, labels=tissue, gene.selection = "common", col=col.tissue, dim.plot = c(1,2))
# title(main="placentas and deciduas
#       common gene selection
#       all transcripts")
# dev.off()
# dev.off()x
# 
# png(filename ="FIGURES/MDS/placentas_deciduas_excludeRemovals_common_top100transcripts_dim1&2.png",width=5,height=5,units="in",res=1200)
# plotMDS(lcpm_tran, labels=tissue, gene.selection = "common", top = 100, col=col.tissue, dim.plot = c(1,2))
# title(main="placentas and deciduas
#       common gene selection
#       top 100 genes")
# dev.off()
# dev.off()
#-----------------------------------------------------

png(filename ="FIGURES/MDS/placentas_deciduas_dim1&2_PAR2x2.png",width=5,height=5,units="in",res=1200)
par(mfrow = c(2,2), mar = c(1,1,1,1) + 0.9)

plotMDS(lcpm_gene, labels=tissue, gene.selection = "common", col=col.tissue, dim.plot = c(1,2))
title(main="All genes")

plotMDS(lcpm_gene, labels=tissue, gene.selection = "common", top = 100, col=col.tissue, dim.plot = c(1,2))
title(main="Top 100 genes")

plotMDS(lcpm_tran, labels=tissue, gene.selection = "common", col=col.tissue, dim.plot = c(1,2))
title(main="All transcripts")

plotMDS(lcpm_tran, labels=tissue, gene.selection = "common", top = 100, col=col.tissue, dim.plot = c(1,2))
title(main="Top 100 transcripts")

dev.off()
dev.off()

