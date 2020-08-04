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
#--------------
# set working directory 
#--------------
setwd("~/Dropbox (ASU)/Placenta/BATCH2_PLACENTA_DECIDUA_ANALYSIS/HISAT_FeatureCounts/batch1_and_batch2/")
#-----------
# Read in count, genes, and phenotype data
#-----------
counts_g <- read.delim("counts_pheno/batch1_and_batch2_geneCounts.tsv", header=TRUE, sep="\t")
colnames(counts_g) <- str_replace_all(colnames(counts_g), pattern="\\.","-") # replace . with - in sample names
genes <- read.csv("counts_pheno/genesID.csv", header=TRUE, sep = "\t")

counts_t <- read.delim("counts_pheno/batch1_and_batch2_transcriptCounts.tsv", header=TRUE, sep="\t")
colnames(counts_t) <- str_replace_all(colnames(counts_t), pattern="\\.","-") # replace . with - in sample names
transcripts <- read.csv("counts_pheno/transcriptsID.csv", header=TRUE, sep = "\t")

pheno <- read.csv("counts_pheno/batch1_and_batch2_FandM_pheno.csv", header=TRUE, sep=",")

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
all_removals <- c(placenta_batch1_removals, placenta_batch2_removals)

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
race<- factor(pheno_ExRemovals$race, levels=c("Asian", "Black", "White", "Other"))
site<- factor(pheno_ExRemovals$site, levels=c("RNA_1", "RNA_2"))
lane<- factor(pheno_ExRemovals$lane, levels=c("L001", "L002","L003", "L004","L005", "L006"))
batch<- factor(pheno_ExRemovals$batch, levels=c("1", "2"))
rep <- factor(pheno_ExRemovals$REP, level=c("OBG0044-1", "OBG0053-1", "OBG0068-1", "OBG0111-1", 
                                            "OBG0112-1", "OBG0115-1", "OBG0116-1", "OBG0117-1", 
                                            "OBG0118-1", "OBG0120-1", "OBG0122-1", "OBG0123-1", 
                                            "OBG0126-1", "OBG0130-1", "OBG0132-1", "OBG0133-1", 
                                            "OBG0156-1", "OBG0158-1", "OBG0166-1", "OBG0170-1", 
                                            "OBG0174-1", "OBG0175-1", "OBG0178-1", "YPOPS0006-1",
                                            "OBG0044", "OBG0053", "OBG0068", "OBG0111", "OBG0112", 
                                            "OBG0115", "OBG0116", "OBG0117", "OBG0118", "OBG0120", 
                                            "OBG0122", "OBG0123", "OBG0126", "OBG0130", "OBG0132", 
                                            "OBG0133", "OBG0156", "OBG0158", "OBG0166", "OBG0170", 
                                            "OBG0174", "OBG0175", "OBG0178", "YPOPS0006", "OBG0014", 
                                            "OBG0015", "OBG0019", "OBG0021", "OBG0022", "OBG0024", 
                                            "OBG0026", "OBG0027", "OBG0028", "OBG0029", "OBG0030", 
                                            "OBG0031", "OBG0032", "OBG0039", "OBG0047", "OBG0050", 
                                            "OBG0051", "OBG0053B2", "OBG0065", "OBG0066", "OBG0085", 
                                            "OBG0090", "OBG0107", "OBG0121", "OBG0138", "OBG0149", 
                                            "OBG0180", "OBG0188", "OBG0191", "OBG0201", "OBG0205", 
                                            "OBG0289", "OBG0338", "OBG0342", "YPOPS0007M", "YPOPS0123M"))

# add groups to samples in dge list 
dge_g$samples$sex <- sex
dge_g$samples$race <- race
dge_g$samples$site <- site
dge_g$samples$lane <- lane
dge_g$samples$batch <- batch
dge_g$samples$rep <- rep
dge_g$genes <- genes


dge_t$samples$sex <- sex
dge_t$samples$race <- race
dge_t$samples$site <- site
dge_t$samples$lane <- lane
dge_t$samples$batch <- batch
dge_t$samples$rep <- rep
dge_t$genes <- genes

#-----------
# data pre-processing
#-----------
dge_g <-sumTechReps(dge_g, dge_g$samples$rep) # comment out to not sum replicates
dge_t <-sumTechReps(dge_t, dge_t$samples$rep) # comment out to not sum replicates

cpm_g <- cpm(dge_g) #log2-transformed counts per gene per million mapped reads (cpm).
lcpm_g <- cpm(dge_g, log=TRUE) # log transformting the data 

cpm_t <- cpm(dge_t) #log2-transformed counts per gene per million mapped reads (cpm).
lcpm_t <- cpm(dge_t, log=TRUE) # log transformting the data 

# filter to only include counts of 1 or greater in at least half of the samples 
keep<- rowSums(cpm_g>1)>=0
table(keep)
dge_g_CPM <- dge_g[keep,, keep.lib.sizes=FALSE]
dim(dge_g_CPM)
logdge_g_CPM1 <- cpm(dge_g_CPM, log=TRUE) # log transformting the data 
dim(logdge_g_CPM1)

cpmBatch_gene <- removeBatchEffect(logdge_g_CPM1, batch=dge_g$samples$batch)

sex_g <- dge_g$samples$sex
batch_g <- dge_g$samples$batch


keep<- rowSums(cpm_t>1)>=0
table(keep)
dge_t_CPM <- dge_t[keep,, keep.lib.sizes=FALSE]
dim(dge_t_CPM)
logdge_t_CPM1 <- cpm(dge_t_CPM, log=TRUE) # log transformting the data 
dim(logdge_t_CPM1)

cpmBatch_tran <- removeBatchEffect(logdge_t_CPM1, batch=dge_t$samples$batch)

sex <- dge_t$samples$sex
batch <- dge_t$samples$batch

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
#---- 
png(filename ="FIGURES/MDS/placentas_batch_1and2/label_batch/placentas_batch_1and2_excludingFailes_common_allgenes_dim1&2.png",width=6,height=6,units="in",res=1200)
plotMDS(logdge_g_CPM1, label=batch, col=col.sex, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="All genes")
dev.off()
dev.off()

png(filename ="FIGURES/MDS/placentas_batch_1and2/label_batch/placentas_batch_1and2_excludingFailes_common_top100genes_dim1&2.png",width=6,height=6,units="in",res=1200)
plotMDS(logdge_g_CPM1, label=batch, col=col.sex, 
        top = 100, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="Top 100 genes")
dev.off()
dev.off()

png(filename ="FIGURES/MDS/placentas_batch_1and2/label_batch/placentas_batch_1and2_excludingFailes_common_allTrans_dim1&2.png",width=6,height=6,units="in",res=1200)
plotMDS(logdge_t_CPM1, label=batch, col=col.sex, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="All transcripts")
dev.off()
dev.off()

png(filename ="FIGURES/MDS/placentas_batch_1and2/label_batch/placentas_batch_1and2_excludingFailes_common_top100trans_dim1&2.png",width=6,height=6,units="in",res=1200)
plotMDS(logdge_t_CPM1, label=batch, col=col.sex, 
        top = 100, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="Top 100 transcripts")
dev.off()
dev.off()


png(filename ="FIGURES/MDS/placentas_batch_1and2/label_batch/placentas_batch_1and2_excludingFails_dim1&2_PAR2x2.png",width=6,height=6,units="in",res=1200)
par(mfrow = c(2,2), mar = c(1,1,1,1) + 0.9)
plotMDS(logdge_g_CPM1, label=batch, col=col.sex, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="All genes")

plotMDS(logdge_g_CPM1, label=batch, col=col.sex, 
        top = 100, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="Top 100 genes")

plotMDS(logdge_t_CPM1, label=batch, col=col.sex, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="All transcripts")

plotMDS(logdge_t_CPM1, label=batch, col=col.sex, 
        top = 100, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="Top 100 transcripts")
dev.off()
dev.off()

#---- 
# label sex, no batch correction
#---- 
png(filename ="FIGURES/MDS/placentas_batch_1and2/label_sex/placentas_batch_1and2_excludingFailes_common_allgenes_dim1&2.png",width=6,height=6,units="in",res=1200)
plotMDS(logdge_g_CPM1, label=sex, col=col.sex, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="All genes")
dev.off()
dev.off()

png(filename ="FIGURES/MDS/placentas_batch_1and2/label_sex/placentas_batch_1and2_excludingFailes_common_top100genes_dim1&2.png",width=6,height=6,units="in",res=1200)
plotMDS(logdge_g_CPM1, label=sex, col=col.sex, 
        top = 100, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="Top 100 genes")
dev.off()
dev.off()

png(filename ="FIGURES/MDS/placentas_batch_1and2/label_sex/placentas_batch_1and2_excludingFailes_common_allTrans_dim1&2.png",width=6,height=6,units="in",res=1200)
plotMDS(logdge_t_CPM1, label=sex, col=col.sex, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="All transcripts")
dev.off()
dev.off()

png(filename ="FIGURES/MDS/placentas_batch_1and2/label_sex/placentas_batch_1and2_excludingFailes_common_top100trans_dim1&2.png",width=6,height=6,units="in",res=1200)
plotMDS(logdge_t_CPM1, label=sex, col=col.sex, 
        top = 100, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="Top 100 transcripts")
dev.off()
dev.off()


png(filename ="FIGURES/MDS/placentas_batch_1and2/label_sex/placentas_batch_1and2_excludingFails_dim1&2_PAR2x2.png",width=6,height=6,units="in",res=1200)
par(mfrow = c(2,2), mar = c(1,1,1,1) + 0.9)
plotMDS(logdge_g_CPM1, label=sex, col=col.sex, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="All genes")

plotMDS(logdge_g_CPM1, label=sex, col=col.sex, 
        top = 100, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="Top 100 genes")

plotMDS(logdge_t_CPM1, label=sex, col=col.sex, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="All transcripts")

plotMDS(logdge_t_CPM1, label=sex, col=col.sex, 
        top = 100, 
        gene.selection = "common", dim.plot = c(1,2))
title(main="Top 100 transcripts")
dev.off()
dev.off()


