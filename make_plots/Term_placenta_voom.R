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
library(tidyverse)
library(ggrepel)
library(variancePartition)
library(doParallel)

#--------------
# set working directory 
#--------------
setwd("~/Dropbox (ASU)/Placenta/BATCH2_PLACENTA_DECIDUA_ANALYSIS/HISAT_FeatureCounts/batch1_and_batch2/")
#-----------
# Read in count, genes, and phenotype data
#-----------
counts_g <- read.delim("counts_pheno/placenta_batch1and2_geneCounts.tsv", header=TRUE, sep="\t")
colnames(counts_g) <- str_replace_all(colnames(counts_g), pattern="\\.","-") # replace . with - in sample names
genes <- read.csv("counts_pheno/genesID.csv", header=TRUE, sep = ",")

counts_t <- read.delim("counts_pheno/placenta_batch1and2_transcriptCounts.tsv", header=TRUE, sep="\t")
colnames(counts_t) <- str_replace_all(colnames(counts_t), pattern="\\.","-") # replace . with - in sample names
transcripts <- read.csv("counts_pheno/transcriptsID.csv", header=TRUE, sep = ",")

#pheno <- read.csv("counts_pheno/placenta_pheno.txt", header=TRUE, sep="\t")
pheno <- read.csv("counts_pheno/200508_placentas_pheno.csv", header=TRUE, sep=",")

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
                              "OBG0026-1", "OBG0026-2", "YPOPS0007M-1","YPOPS0007M-2",
                              "OBG0019-1", "OBG0019-2", "OBG0021-1", "OBG0021-2")

placenta_removals <- c(placenta_batch1_removals, placenta_batch2_removals)

no_removals <- c()
all_removals <- c(placenta_removals)

samplesToRemove <- c(all_removals) # update depending on comparison being made
SAMPLE_LENGTH <- as.numeric(length(samplesToRemove)) # to call later 
half_sample_length <- SAMPLE_LENGTH/2 # half the sample length

removals_g <- (names(counts_g) %in% samplesToRemove[1:SAMPLE_LENGTH]) # for matching names create a value of TRUE or FALSES 
counts_ExRemovals_g <-counts_g[!removals_g] # create a new counts file that excludes (Ex) the removals IDs 

removals_t <- (names(counts_t) %in% samplesToRemove[1:SAMPLE_LENGTH]) # for matching names create a value of TRUE or FALSES 
counts_ExRemovals_t <-counts_t[!removals_t] # create a new counts file that excludes (Ex) the removals IDs 
pheno_ExRemovals <- pheno[! pheno$sample %in% samplesToRemove[1:SAMPLE_LENGTH],] # update 1:16 depending on size of samples to remove

#-----------
# create a DGElist 
#-----------
dge_g <- DGEList(counts=counts_ExRemovals_g, genes=genes)
dge_t <- DGEList(counts=counts_ExRemovals_t, genes=transcripts)
dim(dge_g)
dim(dge_t)
#-----------
# organize sample information 
#-----------
samplenames <- (pheno_ExRemovals$sample) # sample names are not unique because a sample may belong to multiple groups
colnames(dge_g) <- samplenames
colnames(dge_t) <- samplenames

# create groups for the samples  
sex <- factor(pheno_ExRemovals$sex, levels=c("female", "male"))
Reported_race<- factor(pheno_ExRemovals$Reported_race, levels=c("Asian", "Black", "White", "Other"))
site<- factor(pheno_ExRemovals$site, levels=c("RNA_1", "RNA_2"))
Lane<- factor(pheno_ExRemovals$Lane, levels=c("L001", "L002","L003", "L004","L005", "L006"))
batch<- factor(pheno_ExRemovals$batch, levels=c("1", "2"))
rep <- factor(pheno_ExRemovals$REP) 
# level=c("OBG0044-1", "OBG0053-1", "OBG0068-1", "OBG0111-1", 
#                                             "OBG0112-1", "OBG0115-1", "OBG0116-1", "OBG0117-1", 
#                                             "OBG0118-1", "OBG0120-1", "OBG0122-1", "OBG0123-1", 
#                                             "OBG0126-1", "OBG0130-1", "OBG0132-1", "OBG0133-1", 
#                                             "OBG0156", "OBG0158-1", "OBG0166-1", "OBG0170-1", 
#                                             "OBG0174", "OBG0175", "OBG0178", "YPOPS0006",
#                                             "OBG0044", "OBG0053", "OBG0068", "OBG0111", "OBG0112", 
#                                             "OBG0115", "OBG0116", "OBG0117", "OBG0118", "OBG0120", 
#                                             "OBG0122", "OBG0123", "OBG0126", "OBG0130", "OBG0132", 
#                                             "OBG0133", "OBG0156", "OBG0158", "OBG0166", "OBG0170", 
#                                             "OBG0174", "OBG0175", "OBG0178", "YPOPS0006", "OBG0014", 
#                                             "OBG0015", "OBG0019", "OBG0021", "OBG0022", "OBG0024", 
#                                             "OBG0026", "OBG0027", "OBG0028", "OBG0029", "OBG0030", 
#                                             "OBG0031", "OBG0032", "OBG0039", "OBG0047", "OBG0050", 
#                                             "OBG0051", "OBG0053B2", "OBG0065", "OBG0066", "OBG0085", 
#                                             "OBG0090", "OBG0107", "OBG0121", "OBG0138", "OBG0149", 
#                                             "OBG0180", "OBG0188", "OBG0191", "OBG0201", "OBG0205", 
#                                             "OBG0289", "OBG0338", "OBG0342", "YPOPS0123M"))

Maternal_age <- factor(pheno_ExRemovals$Maternal_age)
Maternal_age <- as.numeric(levels(Maternal_age))[Maternal_age]
Gravidity <- factor(pheno_ExRemovals$Gravidity)
Gravidity <- as.numeric(levels(Gravidity))[Gravidity]
Parity <- factor(pheno_ExRemovals$Parity)
Parity <- as.numeric(levels(Parity))[Parity]
BirthWeight <- factor(pheno_ExRemovals$BirthWeight)
BirthWeight <- as.numeric(levels(BirthWeight))[BirthWeight]
Pre.Pregnancy_BMI <- factor(pheno_ExRemovals$Pre.Pregnancy_BMI)
Pre.Pregnancy_BMI <- as.numeric(levels(Pre.Pregnancy_BMI))[Pre.Pregnancy_BMI]
MethodConception <- factor(pheno_ExRemovals$MethodConception)
Ancestry_prediction <- factor(pheno_ExRemovals$Ancestry_prediction)
PC1 <- factor(pheno_ExRemovals$PC1)
PC1 <- as.numeric(levels(PC1))[PC1]
PC2 <- factor(pheno_ExRemovals$PC2)
PC2 <- as.numeric(levels(PC2))[PC2]
Gestational_age <- factor(pheno_ExRemovals$Gestational_age)
Gestational_age <- as.numeric(levels(Gestational_age))[Gestational_age]
# add groups to samples in dge list 
dge_g$samples$sex <- sex
dge_g$samples$Reported_race <- Reported_race
dge_g$samples$site <- site
dge_g$samples$Lane <- Lane
dge_g$samples$batch <- batch
dge_g$samples$rep <- rep
dge_g$samples$Maternal_age <- Maternal_age
dge_g$samples$Gravidity <- Gravidity
dge_g$samples$Parity <- Parity
dge_g$samples$BirthWeight <- BirthWeight
dge_g$samples$Pre.Pregnancy_BMI <- Pre.Pregnancy_BMI
dge_g$samples$MethodConception <- MethodConception
dge_g$samples$Ancestry_prediction <- Ancestry_prediction
dge_g$samples$PC1 <- PC1
dge_g$samples$PC2 <- PC2
dge_g$samples$Gestational_age <- Gestational_age
dge_g$genes <- genes


dge_t$samples$sex <- sex
dge_t$samples$Reported_race <- Reported_race
dge_t$samples$site <- site
dge_t$samples$Lane <- Lane
dge_t$samples$batch <- batch
dge_t$samples$rep <- rep
dge_t$samples$Maternal_age <- Maternal_age
dge_t$samples$Gravidity <- Gravidity
dge_t$samples$Parity <- Parity
dge_t$samples$BirthWeight <- BirthWeight
dge_t$samples$Ancestry_prediction <- Ancestry_prediction
dge_t$samples$Pre.Pregnancy_BMI <- Pre.Pregnancy_BMI
dge_t$samples$MethodConception <- MethodConception
dge_t$samples$PC1 <- PC1
dge_t$samples$PC2 <- PC2
dge_t$samples$Gestational_age <- Gestational_age
dge_t$genes <- transcripts

#-----------
# data pre-processing
#-----------
dge_g <-sumTechReps(dge_g, dge_g$samples$rep) # comment out to not sum replicates
dge_t <-sumTechReps(dge_t, dge_t$samples$rep) # comment out to not sum replicates

#---- 
cpm_gene <- cpm(dge_g) #log2-transformed counts per gene per million mapped reads (cpm).

may <- cbind(genes, cpm_gene)
may$length <- NULL
may2 <- melt(may, id=c("Geneid", "chr"))
pheno_p <- subset(pheno_ExRemovals, site == "RNA_1")
df_merged <- merge(may2, pheno_p, by.x="variable", by.y="REP", all.x=TRUE)
write.table(df_merged, "genelists/placentas_batch_1and2/cpm_placenta_df_genes.txt", sep = "\t")
head(df_merged)

fpkm <- rpkm(dge_g, gene.length=dge_g$genes$Length)
dim(fpkm)
placenta_FPKM <- cbind(genes, fpkm)
write.table(placenta_FPKM, "genelists/placentas_batch_1and2/placenta_geneLevel_FPKM.txt", sep = "\t", quote = FALSE, row.names = FALSE)
female_mean_fpkm <- apply(as.data.frame(fpkm)
                          [(dge_g$samples$sex=="female")],
                          1, mean, na.rm=TRUE)
male_mean_fpkm <- apply(as.data.frame(fpkm)
                        [(dge_g$samples$sex=="male")],
                        1, mean, na.rm=TRUE)
#keep <- (female_mean_fpkm > 0 | male_mean_fpkm > 0)

keep <- (female_mean_fpkm > 1.0 | male_mean_fpkm > 1.0)
#------ 
dge_g <- dge_g[keep,,keep.lib.sizes=FALSE]
dge_g <- calcNormFactors(dge_g, method="TMM")
dim(dge_g$genes)

cpm_gene <- cpm(dge_g) #log2-transformed counts per gene per million mapped reads (cpm).
lcpm_gene <- cpm(cpm_gene, log=TRUE) 
DGE_lcpm <- cbind(dge_g$genes, lcpm_gene)
DGE_cpm <- cbind(dge_g$genes, cpm_gene)

write.table(DGE_cpm, "genelists/placentas_batch_1and2/placenta_cpm_gene.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(DGE_lcpm, "genelists/placentas_batch_1and2/placenta_lcpm_gene.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#---- 
cpm_tran <- cpm(dge_t) #log2-transformed counts per tran per million mapped reads (cpm).

fpkm <- rpkm(dge_t, gene.length=dge_t$genes$Length)
placenta_FPKM <- cbind(transcripts, fpkm)
write.table(placenta_FPKM, "genelists/placentas_batch_1and2/placenta_transLevel_FPKM.txt", sep = "\t", quote = FALSE, row.names = FALSE)

female_mean_fpkm <- apply(as.data.frame(fpkm)
                          [(dge_t$samples$sex=="female")],
                          1, mean, na.rm=TRUE)
male_mean_fpkm <- apply(as.data.frame(fpkm)
                        [(dge_t$samples$sex=="male")],
                        1, mean, na.rm=TRUE)


keep <- (female_mean_fpkm > 1.0 | male_mean_fpkm > 1.0)
#------ 
dge_t <- dge_t[keep,,keep.lib.sizes=FALSE]
dge_t <- calcNormFactors(dge_t, method="TMM")
dim(dge_t$genes)

cpm_tran <- cpm(dge_t) #log2-transformed counts per gene per million mapped reads (cpm).
lcpm_tran <- cpm(cpm_tran, log=TRUE) 

#---------

cpmBatch_gene <- removeBatchEffect(lcpm_gene, batch=dge_g$samples$batch)
cpmBatch_tran <- removeBatchEffect(lcpm_tran, batch=dge_t$samples$batch)

sex <- dge_t$samples$sex
batch <- dge_g$samples$batch
#----------
col.sex <- sex
levels(col.sex) <-  brewer.pal(nlevels(col.sex), "Set2")
col.sex <- as.character(col.sex)

#---- 
# label batch
#---- 
png(filename ="FIGURES/MDS/placentas_batch_1and2/label_batch/placentas_batch_1and2_excludingFailes_batchCORR_common_allgenes_dim1&2.png",width=6,height=6,units="in",res=1200)
plotMDS(cpmBatch_gene, label=batch, col=col.sex, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="All genes")
dev.off()
dev.off()

png(filename ="FIGURES/MDS/placentas_batch_1and2/label_batch/placentas_batch_1and2_excludingFailes_batchCORR_common_top100genes_dim1&2.png",width=6,height=6,units="in",res=1200)
plotMDS(cpmBatch_gene, label=batch, col=col.sex, 
        top = 100, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="Top 100 genes")
dev.off()
dev.off()

png(filename ="FIGURES/MDS/placentas_batch_1and2/label_batch/placentas_batch_1and2_excludingFailes_batchCORR_common_allTrans_dim1&2.png",width=6,height=6,units="in",res=1200)
plotMDS(cpmBatch_tran, label=batch, col=col.sex, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="All transcripts")
dev.off()
dev.off()

png(filename ="FIGURES/MDS/placentas_batch_1and2/label_batch/placentas_batch_1and2_excludingFailes_batchCORR_common_top100trans_dim1&2.png",width=6,height=6,units="in",res=1200)
plotMDS(cpmBatch_tran, label=batch, col=col.sex, 
        top = 100, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="Top 100 transcripts")
dev.off()
dev.off()

png(filename ="FIGURES/MDS/placentas_batch_1and2/label_batch/placentas_batch_1and2_excludingFails_batchCORR_dim1&2_PAR2x2.png",width=6,height=6,units="in",res=1200)
par(mfrow = c(2,2), mar = c(1,1,1,1) + 0.9)
plotMDS(cpmBatch_gene, label=batch, col=col.sex, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="All genes")

plotMDS(cpmBatch_gene, label=batch, col=col.sex, 
        top = 100, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="Top 100 genes")

plotMDS(cpmBatch_tran, label=batch, col=col.sex, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="All transcripts")

plotMDS(cpmBatch_tran, label=batch, col=col.sex, 
        top = 100, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="Top 100 transcripts")
dev.off()
dev.off()

#---- 
# label sex
#---- 
png(filename ="FIGURES/MDS/placentas_batch_1and2/label_sex/placentas_batch_1and2_excludingFailes_batchCORR_common_allgenes_dim1&2.png",width=6,height=6,units="in",res=1200)
plotMDS(cpmBatch_gene, label=sex, col=col.sex, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="All genes")
dev.off()
dev.off()

png(filename ="FIGURES/MDS/placentas_batch_1and2/label_sex/placentas_batch_1and2_excludingFailes_batchCORR_common_top100genes_dim1&2.png",width=6,height=6,units="in",res=1200)
plotMDS(cpmBatch_gene, label=sex, col=col.sex, 
        top = 100, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="Top 100 genes")
dev.off()
dev.off()

png(filename ="FIGURES/MDS/placentas_batch_1and2/label_sex/placentas_batch_1and2_excludingFailes_batchCORR_common_allTrans_dim1&2.png",width=6,height=6,units="in",res=1200)
plotMDS(cpmBatch_tran, label=sex, col=col.sex, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="All transcripts")
dev.off()
dev.off()

png(filename ="FIGURES/MDS/placentas_batch_1and2/label_sex/placentas_batch_1and2_excludingFailes_batchCORR_common_top100trans_dim1&2.png",width=6,height=6,units="in",res=1200)
plotMDS(cpmBatch_tran, label=sex, col=col.sex, 
        top = 100, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="Top 100 transcripts")
dev.off()
dev.off()


png(filename ="FIGURES/MDS/placentas_batch_1and2/label_sex/placentas_batch_1and2_excludingFails_batchCORR_dim1&2_PAR2x2.png",width=6,height=6,units="in",res=1200)
par(mfrow = c(2,2), mar = c(1,1,1,1) + 0.9)
plotMDS(cpmBatch_gene, label=sex, col=col.sex, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="All genes")

plotMDS(cpmBatch_gene, label=sex, col=col.sex, 
        top = 100, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="Top 100 genes")

plotMDS(cpmBatch_tran, label=sex, col=col.sex, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="All transcripts")

plotMDS(cpmBatch_tran, label=sex, col=col.sex, 
        top = 100, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="Top 100 transcripts")
dev.off()
dev.off()

#-----------
# No batch correction
#----------
#---- 
# label batch, no batch correction
# #---- 
# png(filename ="FIGURES/MDS/placentas_batch_1and2/label_batch/placentas_batch_1and2_excludingFailes_common_allgenes_dim1&2.png",width=6,height=6,units="in",res=1200)
# plotMDS(logdge_g_CPM1, label=batch_g, col=col.sex, 
#         gene.selection = "common", dim.plot = c(1,2))
# title(main="All genes")
# dev.off()
# dev.off()
# 
# png(filename ="FIGURES/MDS/placentas_batch_1and2/label_batch/placentas_batch_1and2_excludingFailes_common_top100genes_dim1&2.png",width=6,height=6,units="in",res=1200)
# plotMDS(logdge_g_CPM1, label=batch, col=col.sex, 
#         top = 100, 
#         gene.selection = "common", dim.plot = c(1,2))
# title(main="Top 100 genes")
# dev.off()
# dev.off()
# 
# png(filename ="FIGURES/MDS/placentas_batch_1and2/label_batch/placentas_batch_1and2_excludingFailes_common_allTrans_dim1&2.png",width=6,height=6,units="in",res=1200)
# plotMDS(logdge_t_CPM1, label=batch, col=col.sex, 
#         gene.selection = "common", dim.plot = c(1,2))
# title(main="All transcripts")
# dev.off()
# dev.off()
# 
# png(filename ="FIGURES/MDS/placentas_batch_1and2/label_batch/placentas_batch_1and2_excludingFailes_common_top100trans_dim1&2.png",width=6,height=6,units="in",res=1200)
# plotMDS(logdge_t_CPM1, label=batch, col=col.sex, 
#         top = 100, 
#         gene.selection = "common", dim.plot = c(1,2))
# title(main="Top 100 transcripts")
# dev.off()
# dev.off()
# 
# 
# png(filename ="FIGURES/MDS/placentas_batch_1and2/label_batch/placentas_batch_1and2_excludingFails_dim1&2_PAR2x2.png",width=6,height=6,units="in",res=1200)
# par(mfrow = c(2,2), mar = c(1,1,1,1) + 0.9)
# plotMDS(logdge_g_CPM1, label=batch, col=col.sex, 
#         gene.selection = "common", dim.plot = c(1,2))
# title(main="All genes")
# 
# plotMDS(logdge_g_CPM1, label=batch, col=col.sex, 
#         top = 100, 
#         gene.selection = "common", dim.plot = c(1,2))
# title(main="Top 100 genes")
# 
# plotMDS(logdge_t_CPM1, label=batch, col=col.sex, 
#         gene.selection = "common", dim.plot = c(1,2))
# title(main="All transcripts")
# 
# plotMDS(logdge_t_CPM1, label=batch, col=col.sex, 
#         top = 100, 
#         gene.selection = "common", dim.plot = c(1,2))
# title(main="Top 100 transcripts")
# dev.off()
# dev.off()
# 
# #---- 
# # label sex, no batch correction
# #---- 
# png(filename ="FIGURES/MDS/placentas_batch_1and2/label_sex/placentas_batch_1and2_excludingFailes_common_allgenes_dim1&2.png",width=6,height=6,units="in",res=1200)
# plotMDS(logdge_g_CPM1, label=sex, col=col.sex, 
#         gene.selection = "common", dim.plot = c(1,2))
# title(main="All genes")
# dev.off()
# dev.off()
# 
# png(filename ="FIGURES/MDS/placentas_batch_1and2/label_sex/placentas_batch_1and2_excludingFailes_common_top100genes_dim1&2.png",width=6,height=6,units="in",res=1200)
# plotMDS(logdge_g_CPM1, label=sex, col=col.sex, 
#         top = 100, 
#         gene.selection = "common", dim.plot = c(1,2))
# title(main="Top 100 genes")
# dev.off()
# dev.off()
# 
# png(filename ="FIGURES/MDS/placentas_batch_1and2/label_sex/placentas_batch_1and2_excludingFailes_common_allTrans_dim1&2.png",width=6,height=6,units="in",res=1200)
# plotMDS(logdge_t_CPM1, label=sex, col=col.sex, 
#         gene.selection = "common", dim.plot = c(1,2))
# title(main="All transcripts")
# dev.off()
# dev.off()
# 
# png(filename ="FIGURES/MDS/placentas_batch_1and2/label_sex/placentas_batch_1and2_excludingFailes_common_top100trans_dim1&2.png",width=6,height=6,units="in",res=1200)
# plotMDS(logdge_t_CPM1, label=sex, col=col.sex, 
#         top = 100, 
#         gene.selection = "common", dim.plot = c(1,2))
# title(main="Top 100 transcripts")
# dev.off()
# dev.off()
# 
# 
# png(filename ="FIGURES/MDS/placentas_batch_1and2/label_sex/placentas_batch_1and2_excludingFails_dim1&2_PAR2x2.png",width=6,height=6,units="in",res=1200)
# par(mfrow = c(2,2), mar = c(1,1,1,1) + 0.9)
# plotMDS(logdge_g_CPM1, label=sex, col=col.sex, 
#         gene.selection = "common", dim.plot = c(1,2))
# title(main="All genes")
# 
# plotMDS(logdge_g_CPM1, label=sex, col=col.sex, 
#         top = 100, 
#         gene.selection = "common", dim.plot = c(1,2))
# title(main="Top 100 genes")
# 
# plotMDS(logdge_t_CPM1, label=sex, col=col.sex, 
#         gene.selection = "common", dim.plot = c(1,2))
# title(main="All transcripts")
# 
# plotMDS(logdge_t_CPM1, label=sex, col=col.sex, 
#         top = 100, 
#         gene.selection = "common", dim.plot = c(1,2))
# title(main="Top 100 transcripts")
# dev.off()
# dev.off()
# 
#-----------------------
#
#------------------------
#-------------------------------------------
# Voom
#-------------------------------------------
#--------- gene level Voom 
group <-interaction(dge_g$samples$sex) #, DGE$samples$batch, DGE$samples$Lane, DGE$samples$Reported_race) #, DGE$samples$batch, DGE$samples$Reported_race, DGE$samples$site, DGE$samples$Lane, DGE$samples)

PC1 <- dge_g$samples$PC1
PC2 <- dge_g$samples$PC2
BirthWeight <- dge_g$samples$BirthWeight
Lane <- dge_g$samples$Lane
batch <- dge_g$samples$batch

mm <- model.matrix(~0 + group + batch + Lane + PC1 + PC2)
colnames(mm) <- gsub("v\\$targets\\$sex", "", colnames(mm))
head(mm)

y <- voom(dge_g, mm, plot = T) # plot T for plotting 

# Fitting linear models in limma
fit <- lmFit(y, mm)


png(filename ="FIGURES/voomFits/placentas_batch_1and2/batch_lane_PC1_PC2/placenta_voomFit_geneLevel.png",width=6,height=6,units="in",res=1200)
v <- voomWithQualityWeights(dge_g, mm, plot=TRUE)
dev.off()
dev.off()

voomCounts <- v$E
# 2.7325346
voomCountsMatrix <- data.matrix(voomCounts, rownames.force = NA)
genes <- v$genes
fit <- lmFit(v, mm)
fit
contrasts <- makeContrasts(FvsM = groupfemale - groupmale, 
                           levels=colnames(mm))


#-------------------------------
# # variance Part
# # Fit model and extract results. A linear mixed model is used
# # each entry in results is a regression model fit on a single gene
# 
# #form <- ~ Maternal_age + Gestational_age + Gravidity + Parity + Pre.Pregnancy_BMI + BirthWeight + (1|Reported_race) + (1|MethodConception) + (1|Lane) + (1|sex)
# form <- ~ batch + Maternal_age + Gestational_age + Gravidity + Parity + Pre.Pregnancy_BMI + BirthWeight + (1|Reported_race) + (1|MethodConception) + (1|Lane) + (1|sex)
# 
# # 2) extract variance fractions from each model fit for each gene, returns fraction of variation attributable to each variable
# # Interpretation: the variance explained by each variables after correcting for all other variables
# 
# cl <- makeCluster(4)
# registerDoParallel(cl)
# PHENOvarpart <- subset(pheno_ExRemovals, pheno_ExRemovals$site == "RNA_1")
# # Fitting the linear model
# varPart <- fitExtractVarPartModel(voomCounts, form, PHENOvarpart)
# # Error Some predictor variables are on very different scales: consider rescaling 
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
# write.table(genesVP, "variancePart/placentas/varPart_placentas.txt")
# 
# # violin plot of contribution of each variable to total variance
# vp <- sortCols( varPart )
# png(filename ="variancePart/placentas/varPart_placentas.png",width=6,height=6,units="in",res=1200)
# plotVarPart( vp )
# dev.off()
# dev.off()

#---- end var part

# Running contrast analysis
vfit <- contrasts.fit(fit, contrasts = contrasts)
# Looking at N of DEGs with adj. p <0.01 and log2FC>2
sumTable <- summary(decideTests(vfit, adjust.method = "BH", p.value = 0.05, lfc = 1))

# Computing differential expression based on the empirical Bayes moderation of
# the standard errors towards a common value. Robust = should the estimation of
# the empirical Bayes prior parameters be robustified against outlier sample
# variances?
veBayesFit <- eBayes(vfit, robust=TRUE)
png(filename ="FIGURES/voomFits/placentas_batch_1and2/batch_lane_PC1_PC2/placenta_voomFit_geneLevel_FinalModel.png",width=6,height=6,units="in",res=1200)
plotSA(veBayesFit, main = "Final model: Mean-variance trend")
dev.off()
dev.off()

# Getting the DEGs. Log2FC of 1 is equivalent to linear fold change of 2.
# Getting summary statistics for all genes
allComparisons <- colnames(contrasts)
coef = 1

for (i in allComparisons){
  vTopTableAll <- topTable(veBayesFit, coef=coef, n=Inf, p.value=1, lfc=0)
  path <- paste("./DEGs/placentas_batch_1and2/batch_lane_PC1_PC2/DEGs_placentas_genes", i, "_fdr1_lfc0.txt", sep = "")
  write.table(vTopTableAll, path, sep = "\t")
  
  # Adj.p<0.05, log2FC>0
  vTopTableFDR.05 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.05, lfc=0)
  path <- paste("./DEGs/placentas_batch_1and2/batch_lane_PC1_PC2/DEGs_placentas_genes", i, "_fdr05_lfc0.txt", sep = "")
  write.table(vTopTableFDR.05, path, sep = "\t")
  
  # log2FC>1
  vTopTableFC1 <- topTable(veBayesFit, coef=coef, n=Inf, lfc=1)
  path <- paste("./DEGs/placentas_batch_1and2/batch_lane_PC1_PC2/DEGs_placentas_genes", i, "_lfc1.txt", sep = "")
  write.table(vTopTableFC1, path, sep = "\t")
  
  # log2FC>2
  vTopTableFC2 <- topTable(veBayesFit, coef=coef, n=Inf, lfc=2)
  path <- paste("./DEGs/placentas_batch_1and2/batch_lane_PC1_PC2/DEGs_placentas_genes", i, "_lfc2.txt", sep = "")
  write.table(vTopTableFC2, path, sep = "\t")
  
  # Adj.p<0.05, log2FC>1
  vTopTableFDR.05FC1 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.05, lfc=1)
  path <- paste("./DEGs/placentas_batch_1and2/batch_lane_PC1_PC2/DEGs_placentas_genes", i, "_fdr05_lfc1.txt", sep = "")
  write.table(vTopTableFDR.05FC1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>0
  vTopTableFDR.01 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=0)
  path <- paste("./DEGs/placentas_batch_1and2/batch_lane_PC1_PC2/DEGs_placentas_genes", i, "_fdr001_lfc0.txt", sep = "")
  write.table(vTopTableFDR.01, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>1
  vTopTableFDR.01FC1 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=1)
  path <- paste("./DEGs/placentas_batch_1and2/batch_lane_PC1_PC2/DEGs_placentas_genes", i, "_fdr001_lfc1.txt", sep = "")
  write.table(vTopTableFDR.01FC1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>2
  vTopTableFDR.01FC2 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=2)
  path <- paste("./DEGs/placentas_batch_1and2/batch_lane_PC1_PC2/DEGs_placentas_genes", i, "_fdr001_lfc2.txt", sep = "")
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
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") 
#                       | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") 
#                       | adj.P.Val<0.05 & logFC < -1 & (chr=="chrX") & (color=="4") 
#                       |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))

#subset_data <- subset(df.t, (logFC < -1.5 & P.Value <0.01 & color=="1") | (logFC > 1.5 & P.Value <0.01 & color=="1"))

subset_data <- subset(df.t, logFC < -1 & (chr=="chrY") & (color=="5") 
                      | logFC > 1 & (chr=="chrX") & (color=="4") 
                      | logFC < -1 & (chr=="chrX") & (color=="4") 
                      | (chr!="chrX") & (chr!="chrY") & (color=="3") & logFC >1
                      | (chr!="chrX") & (chr!="chrY") & (color=="3") & logFC < -1 )
subset_data<- unique(subset_data)

# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/placentas_batch_1and2/batch_lane_PC1_PC2/placenta_MvsF_geneLevel_fdr05_lfc1_adjpsizeSet.png",  width= 440, height= 540)
p2 <- p + ggtitle("")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1, size = 5) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
p2
dev.off()
dev.off()dev.off()
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
png(filename ="Figures/Volcano/placentas_batch_1and2/batch_lane_PC1_PC2/placenta_MvsF_geneLevel_fdr05_lfc1.png", width= 480, height= 580)
p + ggtitle("placenta differential expression gene level")+
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
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/placentas_batch_1and2/batch_lane_PC1_PC2/placenta_MvsF_geneLevel_lfc1_sizeSet.png", width= 480, height= 580)
p + ggtitle("placenta differential expression gene level")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
p
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
png(filename ="Figures/Volcano/placentas_batch_1and2/batch_lane_PC1_PC2/placenta_MvsF_geneLevel_lfc1.png", width= 480, height= 580)
p + ggtitle("placenta differential expression gene level")+
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
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/placentas_batch_1and2/batch_lane_PC1_PC2/placenta_MvsF_geneLevel_fdr05_sizeSet.png", width= 480, height= 580)
p + ggtitle("placenta differential expression gene level")+
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
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/placentas_batch_1and2/batch_lane_PC1_PC2/placenta_MvsF_geneLevel_fdr05.png", width= 480, height= 580)
p + ggtitle("placenta differential expression gene level")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()


#--------- transcript level Voom 

#--------- transcript level Voom 
group <-interaction(dge_t$samples$sex) #, DGE$samples$batch, DGE$samples$Lane, DGE$samples$Reported_race) #, DGE$samples$batch, DGE$samples$Reported_race, DGE$samples$site, DGE$samples$Lane, DGE$samples)

PC1 <- dge_t$samples$PC1
PC2 <- dge_t$samples$PC2
BirthWeight <- dge_t$samples$BirthWeight
Lane <- dge_t$samples$Lane
batch <- dge_t$samples$batch

mm <- model.matrix(~0 + group + batch + Lane + PC1 + PC2)
colnames(mm) <- gsub("v\\$targets\\$sex", "", colnames(mm))
head(mm)

y <- voom(dge_t, mm, plot = T) # plot T for plotting 

# Fitting linear models in limma
fit <- lmFit(y, mm)

png(filename ="FIGURES/voomFits/placentas_batch_1and2/batch_lane_PC1_PC2/placenta_voomFit_transLevel.png",width=6,height=6,units="in",res=1200)
v <- voomWithQualityWeights(dge_t, mm, plot=TRUE)
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
png(filename ="FIGURES/voomFits/placentas_batch_1and2/batch_lane_PC1_PC2/placenta_voomFit_transLevel_FinalModel.png",width=6,height=6,units="in",res=1200)
plotSA(veBayesFit, main = "Final model: Mean-variance trend")
dev.off()
dev.off()

# Getting the DEGs. Log2FC of 1 is equivalent to linear fold change of 2.
# Getting summary statistics for all genes
allComparisons <- colnames(contrasts)
coef = 1

for (i in allComparisons){
  vTopTableAll <- topTable(veBayesFit, coef=coef, n=Inf, p.value=1, lfc=0)
  path <- paste("./DEGs/placentas_batch_1and2/batch_lane_PC1_PC2/DEGs_placentas_trans", i, "_fdr1_lfc0.txt", sep = "")
  write.table(vTopTableAll, path, sep = "\t")
  
  # Adj.p<0.05, log2FC>0
  vTopTableFDR.05 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.05, lfc=0)
  path <- paste("./DEGs/placentas_batch_1and2/batch_lane_PC1_PC2/DEGs_placentas_trans", i, "_fdr05_lfc0.txt", sep = "")
  write.table(vTopTableFDR.05, path, sep = "\t")
  
  # log2FC>1
  vTopTableFC1 <- topTable(veBayesFit, coef=coef, n=Inf, lfc=1)
  path <- paste("./DEGs/placentas_batch_1and2/batch_lane_PC1_PC2/DEGs_placentas_trans", i, "_lfc1.txt", sep = "")
  write.table(vTopTableFC1, path, sep = "\t")
  
  # log2FC>2
  vTopTableFC2 <- topTable(veBayesFit, coef=coef, n=Inf, lfc=2)
  path <- paste("./DEGs/placentas_batch_1and2/batch_lane_PC1_PC2/DEGs_placentas_trans", i, "_lfc2.txt", sep = "")
  write.table(vTopTableFC2, path, sep = "\t")
  
  # Adj.p<0.05, log2FC>1
  vTopTableFDR.05FC1 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.05, lfc=1)
  path <- paste("./DEGs/placentas_batch_1and2/batch_lane_PC1_PC2/DEGs_placentas_trans", i, "_fdr05_lfc1.txt", sep = "")
  write.table(vTopTableFDR.05FC1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>0
  vTopTableFDR.01 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=0)
  path <- paste("./DEGs/placentas_batch_1and2/batch_lane_PC1_PC2/DEGs_placentas_trans", i, "_fdr001_lfc0.txt", sep = "")
  write.table(vTopTableFDR.01, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>1
  vTopTableFDR.01FC1 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=1)
  path <- paste("./DEGs/placentas_batch_1and2/batch_lane_PC1_PC2/DEGs_placentas_trans", i, "_fdr001_lfc1.txt", sep = "")
  write.table(vTopTableFDR.01FC1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>2
  vTopTableFDR.01FC2 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=2)
  path <- paste("./DEGs/placentas_batch_1and2/batch_lane_PC1_PC2/DEGs_placentas_trans", i, "_fdr001_lfc2.txt", sep = "")
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
png(filename ="Figures/Volcano/placentas_batch_1and2/batch_lane_PC1_PC2/placenta_MvsF_tranLevel_fdr05_lfc1_adjpsizeSet.png", width= 480, height= 580)
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
png(filename ="Figures/Volcano/placentas_batch_1and2/batch_lane_PC1_PC2/placenta_MvsF_tranLevel_fdr05_lfc1.png", width= 480, height= 580)
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
png(filename ="Figures/Volcano/placentas_batch_1and2/batch_lane_PC1_PC2/placenta_MvsF_tranLevel_lfc1_sizeSet.png", width= 480, height= 580)
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
png(filename ="Figures/Volcano/placentas_batch_1and2/batch_lane_PC1_PC2/placenta_MvsF_tranLevel_lfc1.png", width= 480, height= 580)
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
png(filename ="Figures/Volcano/placentas_batch_1and2/batch_lane_PC1_PC2/placenta_MvsF_tranLevel_fdr05_sizeSet.png", width= 480, height= 580)
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
png(filename ="Figures/Volcano/placentas_batch_1and2/batch_lane_PC1_PC2/placenta_MvsF_tranLevel_fdr05.png", width= 480, height= 580)
p + ggtitle("placenta differential expression gene level")+
  #geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene_name), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()

#-------------------------------------------
# Voom
#-------------------------------------------
#--------- gene level Voom 
group <-interaction(dge_g$samples$sex) #, DGE$samples$batch, DGE$samples$Lane, DGE$samples$Reported_race) #, DGE$samples$batch, DGE$samples$Reported_race, DGE$samples$site, DGE$samples$Lane, DGE$samples)

PC1 <- dge_g$samples$PC1
PC2 <- dge_g$samples$PC2
BirthWeight <- dge_g$samples$BirthWeight
Lane <- dge_g$samples$Lane
batch <- dge_g$samples$batch

mm <- model.matrix(~0 + group + batch + BirthWeight + Lane + PC1 + PC2)
colnames(mm) <- gsub("v\\$targets\\$sex", "", colnames(mm))
head(mm)

y <- voom(dge_g, mm, plot = T) # plot T for plotting 

# Fitting linear models in limma
fit <- lmFit(y, mm)

png(filename ="FIGURES/voomFits/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/placenta_voomFit_geneLevel.png",width=6,height=6,units="in",res=1200)
v <- voomWithQualityWeights(dge_g, mm, plot=TRUE)
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
pheno_p <- subset(pheno_ExRemovals, site == "RNA_1")
df_merged <- merge(may2, pheno_p, by.x="variable", by.y="REP", all.x=TRUE)
write.table(df_merged, "genelists/placentas_batch_1and2/vfit_placenta_df_genes.txt", sep = "\t")


# Computing differential expression based on the empirical Bayes moderation of
# the standard errors towards a common value. Robust = should the estimation of
# the empirical Bayes prior parameters be robustified against outlier sample
# variances?
veBayesFit <- eBayes(vfit, robust=TRUE)
png(filename ="FIGURES/voomFits/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/placenta_voomFit_geneLevel_FinalModel.png",width=6,height=6,units="in",res=1200)
plotSA(veBayesFit, main = "Final model: Mean-variance trend")
dev.off()
dev.off()

# Getting the DEGs. Log2FC of 1 is equivalent to linear fold change of 2.
# Getting summary statistics for all genes
allComparisons <- colnames(contrasts)
coef = 1

# vTopTableAll <- topTable(veBayesFit, coef=1, n=Inf, p.value=1, lfc=0)
# path <- paste("./DEGs/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/DEGs_placentas_genes_fdr1_lfc0_NOFPKMFILTER.txt", sep = "")
# write.table(vTopTableAll, path, sep = "\t")


for (i in allComparisons){
  vTopTableAll <- topTable(veBayesFit, coef=coef, n=Inf, p.value=1, lfc=0)
  path <- paste("./DEGs/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/DEGs_placentas_genes", i, "_fdr1_lfc0.txt", sep = "")
  write.table(vTopTableAll, path, sep = "\t")
  
  # Adj.p<0.05, log2FC>0
  vTopTableFDR.05 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.05, lfc=0)
  path <- paste("./DEGs/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/DEGs_placentas_genes", i, "_fdr05_lfc0.txt", sep = "")
  write.table(vTopTableFDR.05, path, sep = "\t")
  
  # log2FC>1
  vTopTableFC1 <- topTable(veBayesFit, coef=coef, n=Inf, lfc=1)
  path <- paste("./DEGs/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/DEGs_placentas_genes", i, "_lfc1.txt", sep = "")
  write.table(vTopTableFC1, path, sep = "\t")
  
  # log2FC>2
  vTopTableFC2 <- topTable(veBayesFit, coef=coef, n=Inf, lfc=2)
  path <- paste("./DEGs/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/DEGs_placentas_genes", i, "_lfc2.txt", sep = "")
  write.table(vTopTableFC2, path, sep = "\t")
  
  # Adj.p<0.05, log2FC>1
  vTopTableFDR.05FC1 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.05, lfc=1)
  path <- paste("./DEGs/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/DEGs_placentas_genes", i, "_fdr05_lfc1.txt", sep = "")
  write.table(vTopTableFDR.05FC1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>0
  vTopTableFDR.01 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=0)
  path <- paste("./DEGs/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/DEGs_placentas_genes", i, "_fdr001_lfc0.txt", sep = "")
  write.table(vTopTableFDR.01, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>1
  vTopTableFDR.01FC1 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=1)
  path <- paste("./DEGs/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/DEGs_placentas_genes", i, "_fdr001_lfc1.txt", sep = "")
  write.table(vTopTableFDR.01FC1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>2
  vTopTableFDR.01FC2 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=2)
  path <- paste("./DEGs/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/DEGs_placentas_genes", i, "_fdr001_lfc2.txt", sep = "")
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
  scale_color_manual(values = c("azure3", "azure3", "#66A61E", "#D95F02", "#7570B3")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=18), 
        axis.text.x=element_text(size=18)) +
  theme(axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=18))
p
# Getting a subset of genes to label
# subset_data <- subset(df.t, adj.P.Val<0.05 & logFC < -1 & (chr=="chrY") & (color=="5") | adj.P.Val<0.05 & logFC > 1 & (chr=="chrX") & (color=="4") )
subset_data <- subset(df.t, logFC < -1 & (chr=="chrY") & (color=="5") 
                      | logFC > 1 & (chr=="chrX") & (color=="4") 
                      | logFC < -1 & (chr=="chrX") & (color=="4") 
                      | (chr!="chrX") & (chr!="chrY") & (color=="3") & logFC >1
                      | (chr!="chrX") & (chr!="chrY") & (color=="3") & logFC < -1 )
subset_data<- unique(subset_data)

# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/placenta_MvsF_geneLevel_fdr05_lfc1_adjpsizeSet.png", width= 440, height= 540)
p2 <- p + ggtitle("")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), 
                  segment.alpha = 1, size = 5.5) + 
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
png(filename ="Figures/Volcano/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/placenta_MvsF_geneLevel_fdr05_lfc1.png", width= 480, height= 580)
p + ggtitle("placenta differential expression gene level")+
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
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/placenta_MvsF_geneLevel_lfc1_sizeSet.png", width= 480, height= 580)
p + ggtitle("placenta differential expression gene level")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
p
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
png(filename ="Figures/Volcano/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/placenta_MvsF_geneLevel_lfc1.png", width= 480, height= 580)
p + ggtitle("placenta differential expression gene level")+
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
png(filename ="Figures/Volcano/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/placenta_MvsF_geneLevel_fdr05_sizeSet.png", width= 440, height= 540)
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
subset_data <- subset(df.t, adj.P.Val<0.05  & (chr=="chrY") & (color=="5") 
                      | adj.P.Val<0.05  & (chr=="chrX") & (color=="4") 
                      | adj.P.Val<0.05  & (chr=="chrX") & (color=="4") 
                      |  adj.P.Val<0.05 & (chr!="chrX") & (chr!="chrY") & (color=="3"))
# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
hadjpval <- (-log10(max(df.t$P.Value[df.t$adj.P.Val<0.05], na.rm=TRUE)))
png(filename ="Figures/Volcano/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/placenta_MvsF_geneLevel_fdr05.png", width= 480, height= 580)
p + ggtitle("placenta differential expression gene level")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()


#--------- transcript level Voom 

#--------- transcript level Voom 
group <-interaction(dge_t$samples$sex) #, DGE$samples$batch, DGE$samples$Lane, DGE$samples$Reported_race) #, DGE$samples$batch, DGE$samples$Reported_race, DGE$samples$site, DGE$samples$Lane, DGE$samples)

PC1 <- dge_t$samples$PC1
PC2 <- dge_t$samples$PC2
BirthWeight <- dge_t$samples$BirthWeight
Lane <- dge_t$samples$Lane
batch <- dge_t$samples$batch

mm <- model.matrix(~0 + group + batch + BirthWeight + Lane + PC1 + PC2)
colnames(mm) <- gsub("v\\$targets\\$sex", "", colnames(mm))
head(mm)

y <- voom(dge_t, mm, plot = T) # plot T for plotting 

# Fitting linear models in limma
fit <- lmFit(y, mm)

png(filename ="FIGURES/voomFits/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/placenta_voomFit_transLevel.png",width=6,height=6,units="in",res=1200)
v <- voomWithQualityWeights(dge_t, mm, plot=TRUE)
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
png(filename ="FIGURES/voomFits/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/placenta_voomFit_transLevel_FinalModel.png",width=6,height=6,units="in",res=1200)
plotSA(veBayesFit, main = "Final model: Mean-variance trend")
dev.off()
dev.off()

# Getting the DEGs. Log2FC of 1 is equivalent to linear fold change of 2.
# Getting summary statistics for all genes
allComparisons <- colnames(contrasts)
coef = 1

for (i in allComparisons){
  vTopTableAll <- topTable(veBayesFit, coef=coef, n=Inf, p.value=1, lfc=0)
  path <- paste("./DEGs/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/DEGs_placentas_trans", i, "_fdr1_lfc0.txt", sep = "")
  write.table(vTopTableAll, path, sep = "\t")
  
  # Adj.p<0.05, log2FC>0
  vTopTableFDR.05 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.05, lfc=0)
  path <- paste("./DEGs/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/DEGs_placentas_trans", i, "_fdr05_lfc0.txt", sep = "")
  write.table(vTopTableFDR.05, path, sep = "\t")
  
  # log2FC>1
  vTopTableFC1 <- topTable(veBayesFit, coef=coef, n=Inf, lfc=1)
  path <- paste("./DEGs/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/DEGs_placentas_trans", i, "_lfc1.txt", sep = "")
  write.table(vTopTableFC1, path, sep = "\t")
  
  # log2FC>2
  vTopTableFC2 <- topTable(veBayesFit, coef=coef, n=Inf, lfc=2)
  path <- paste("./DEGs/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/DEGs_placentas_trans", i, "_lfc2.txt", sep = "")
  write.table(vTopTableFC2, path, sep = "\t")
  
  # Adj.p<0.05, log2FC>1
  vTopTableFDR.05FC1 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.05, lfc=1)
  path <- paste("./DEGs/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/DEGs_placentas_trans", i, "_fdr05_lfc1.txt", sep = "")
  write.table(vTopTableFDR.05FC1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>0
  vTopTableFDR.01 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=0)
  path <- paste("./DEGs/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/DEGs_placentas_trans", i, "_fdr001_lfc0.txt", sep = "")
  write.table(vTopTableFDR.01, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>1
  vTopTableFDR.01FC1 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=1)
  path <- paste("./DEGs/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/DEGs_placentas_trans", i, "_fdr001_lfc1.txt", sep = "")
  write.table(vTopTableFDR.01FC1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>2
  vTopTableFDR.01FC2 <- topTable(veBayesFit, coef=coef, n=Inf, p.value=0.01, lfc=2)
  path <- paste("./DEGs/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/DEGs_placentas_trans", i, "_fdr001_lfc2.txt", sep = "")
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
png(filename ="Figures/Volcano/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/placenta_MvsF_tranLevel_fdr05_lfc1_adjpsizeSet.png", width= 480, height= 580)
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
png(filename ="Figures/Volcano/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/placenta_MvsF_tranLevel_fdr05_lfc1.png", width= 480, height= 580)
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
png(filename ="Figures/Volcano/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/placenta_MvsF_tranLevel_lfc1_sizeSet.png", width= 480, height= 580)
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
png(filename ="Figures/Volcano/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/placenta_MvsF_tranLevel_lfc1.png", width= 480, height= 580)
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
png(filename ="Figures/Volcano/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/placenta_MvsF_tranLevel_fdr05_sizeSet.png", width= 480, height= 580)
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
png(filename ="Figures/Volcano/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/placenta_MvsF_tranLevel_fdr05.png", width= 480, height= 580)
p + ggtitle("placenta differential expression gene level")+
  #geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene_name), segment.alpha = 1) + 
  geom_hline(yintercept = hadjpval, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = 1, colour="#A6761D", linetype="dashed"
  ) + geom_vline(xintercept = -1, colour="#A6761D", linetype="dashed")
dev.off()
dev.off()


