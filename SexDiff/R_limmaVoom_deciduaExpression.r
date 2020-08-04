#------------------------------------------------------------------------------------
# Single factor DE analysis using limma/voom
#
#       created by Kimberly Olney, July 19th 2019
#
#------------------------------------------------------------------------------------
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
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")
#install.packages("ggpubr")

library(ggpubr)
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
counts <- read.delim("counts_pheno/placenta_batch2_RNA_HISAT_FC_strandedRF_transcriptCounts_Deciduas.tsv", header = TRUE, sep="\t")
colnames(counts) <- str_replace_all(colnames(counts), pattern="\\.","-") # replace . with - in sample names
#pheno <- read.csv("batch1_and_batch2_FandM_pheno.csv", header=TRUE, sep=",")
pheno <- read.csv("counts_pheno/batch2_deciduas_pheno.csv", header=TRUE, sep=",")
genes <- read.csv("counts_pheno/transcripts_gene_names.csv", header=TRUE, sep = "\t")

decidua_removals <- c("YPOPS0007M-DEC", "OBG0021-DEC")
samplesToRemove <- c(decidua_removals) # update depending on comparison being made
SAMPLE_LENGTH<-as.numeric(length(samplesToRemove)) # to call later 
half_sample_length <- SAMPLE_LENGTH/2 # half the sample length
# update counts data frame to exlude select samples 
removals<-(names(counts) %in% samplesToRemove[1:SAMPLE_LENGTH]) # for matching names create a value of TRUE or FALSES 
counts_ExRemovals <-counts[!removals] # create a new counts file that excludes (Ex) the removals IDs 
geneCounts <- cbind(genes, counts_ExRemovals)

# update pheno data frame to exlude select samples 
pheno_ExRemovals <- pheno[! pheno$sample %in% samplesToRemove[1:SAMPLE_LENGTH],] # update 1:16 depending on size of samples to remove

# create DGE list object 
DGE <- DGEList(counts=counts_ExRemovals, genes = genes) # create DGEList object 
samplenames <- (pheno_ExRemovals$sample) # sample names are not unique because a sample may belong to multiple groups
colnames(DGE) <- samplenames

# create groups for the samples  
#sex <- factor(pheno_ExRemovals$sex, levels=c("female", "male"))
offspring_sex <- factor(pheno_ExRemovals$offspring_sex, levels=c("female", "male"))
#batch <- factor(pheno_ExRemovals$batch, levels=c("1","2"))
race <- factor(pheno_ExRemovals$race, levels=c("Asian", "Black", "White", "Hispanic", "Other"))
#site <- factor(pheno_ExRemovals$site, levels=c("RNA_1", "RNA_2"))
lane <- factor(pheno_ExRemovals$lane, levels=c("L001", "L002","L003", "L004","L005", "L006"))
#rep <- factor(pheno_ExRemovals$REP, level=c("OBG0044", "OBG0053", "OBG0068", "OBG0111", "OBG0112", "OBG0115", "OBG0116", "OBG0117", "OBG0118", "OBG0120", "OBG0122", "OBG0123", "OBG0126", "OBG0130", "OBG0132", "OBG0133", "OBG0156", "OBG0158", "OBG0166", "OBG0170", "OBG0174", "OBG0175", "OBG0178", "YPOPS0006", "OBG0014", "OBG0015", "OBG0019", "OBG0021", "OBG0022", "OBG0024", "OBG0026", "OBG0027", "OBG0028", "OBG0029", "OBG0030", "OBG0031", "OBG0032", "OBG0039", "OBG0047", "OBG0050", "OBG0051", "OBG0053B2", "OBG0065", "OBG0066", "OBG0085", "OBG0090", "OBG0107", "OBG0121", "OBG0138", "OBG0149", "OBG0180", "OBG0188", "OBG0191", "OBG0201", "OBG0205", "OBG0289", "OBG0338", "OBG0342", "YPOPS0007M", "YPOPS0123M"))

# add groups to samples in dge list 
#DGE$samples$sex <- sex
DGE$samples$offspring_sex <- offspring_sex
#DGE$samples$batch <- batch
DGE$samples$race <- race
#DGE$samples$site <- site
DGE$samples$lane <- lane
#DGE$samples$rep <- rep
DGE$genes <- genes

# sum replciates 
#DGE <-sumTechReps(DGE, DGE$samples$rep) # comment out to not sum replicates 
sample_length_sumTech<-length(DGE$samples$group)
half_sample_length_sumTech <- sample_length_sumTech/2
DGE <- calcNormFactors(DGE) # Calculate normalization factors 
# NOTE calcNormFactors doesn’t normalize the data, it just calculates normalization factors for use downstream.
dim(DGE) 

#---- 
cpm <- cpm(DGE) #log2-transformed counts per gene per million mapped reads (cpm).
genesNames <- genes[4]
decidua_cpm <- cbind(genesNames, cpm)
write.table(decidua_cpm, "decidua_cpm.txt", sep ="\t", row.names = FALSE)
keep<- rowSums(cpm>1)>=half_sample_length_sumTech # in at least half the samples 
table(keep)
dge_CPM <- DGE[keep,, keep.lib.sizes=FALSE] # When you subset a DGEList and specify keep.lib.sizes=FALSE, the lib.size for each sample (cf. the y$samples data.frame) will be recalculated to be the sum of the counts left in the rows of the experiment for each sample.
dim(dge_CPM)
logdge_CPM <- cpm(dge_CPM, log=TRUE) # log transformting the data 
dim(logdge_CPM)

# Create a new variable “group” that combines pheno information
#group <-interaction(dge_CPM$samples$sex) #, DGE$samples$batch, DGE$samples$lane, DGE$samples$race) #, DGE$samples$batch, DGE$samples$race, DGE$samples$site, DGE$samples$lane, DGE$samples)
group <-interaction(dge_CPM$samples$offspring_sex) #, DGE$samples$batch, DGE$samples$lane, DGE$samples$race) #, DGE$samples$batch, DGE$samples$race, DGE$samples$site, DGE$samples$lane, DGE$samples)

#batch <- dge_CPM$samples$batch

# Below specifies a model where each coefficient corresponds to a group mean
# Specify the model to be fitted. We do this before using voom since voom uses variances of the model residuals (observed - fitted)
# mm <- model.matrix(~0 + group) # To adjust for batch in the analysis, add batch to the end of the call to model matrix. Everything else about the code stays the same
#mm <- model.matrix(~0 + group + batch) # for including batch in the model 
mm <- model.matrix(~0 + group)
y <- voom(dge_CPM, mm, plot = T) # plot T for plotting 

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
# write.table(FvsM_sig, file="genelists/FvsM_sig_decidua.txt",sep="\t", row.names=FALSE)
write.table(FvsM_sig, file="genelists/FvsM_sig_transcript_decidua.txt",sep="\t", row.names=FALSE)


col.sex <- offspring_sex
levels(col.sex) <-  brewer.pal(nlevels(col.sex), "Set1")
col.sex <- as.character(col.sex)
png(filename ="figures/MDS_decidua_wholeTranscriptome_samples_cpm1_dim1&2_Top20_decidua.png",width=6,height=6,units="in",res=1200)
plotMDS(DGE, labels = offspring_sex, top = 100, gene.selection = "common", col=col.sex, dim.plot = c(1,2))
title(main="MDS decidua by sex of offspring
      whole transcriptome, Top 100")
dev.off()
dev.off()

glMDPlot(tmp, coef=1, status=dt, main=colnames(tmp)[1], counts=dge_CPM$counts, groups=group, launch=TRUE)

# Assigning color to highlight genes that have an absolute fold change > 2 
# and a p-value < Bonferroni cut-off
MvsF_All <- topTable(tmp, n=Inf, coef=1)
P.Value <- c(MvsF_All$P.Value)
logFC <- c(MvsF_All$logFC)
adj.P.Val <- c(MvsF_All$adj.P.Val)
df <- data.frame(P.Value, logFC, adj.P.Val, MvsF_All$gene_name, MvsF_All$chr)

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
colnames(df.t)[4] <- "Geneid"
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
png(filename ="figures/MvsF_VOLCANOplotFC_transcript_decidua.png", width= 480, height= 580)
p + ggtitle("Decidua differential expression")+
  geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=Geneid))+geom_hline(yintercept = 2, colour="red", linetype="dashed"
) + geom_vline(xintercept = 1, colour="red", linetype="dashed"
) + geom_vline(xintercept = -1, colour="red", linetype="dashed")
dev.off()
dev.off()

decF<-subset(pheno_ExRemovals, pheno_ExRemovals$offspring_sex == "female")

# 
#plotMD(tmp, column=1, status=dt[,1], main=colnames(tmp)[1])#, xlim=c(-8,13))
# 
FvsM_topTable <- top.table$Geneid[1:20]
i <- which(y$genes$Geneid %in% FvsM_topTable)
mycol <- colorpanel(1000,"green","white","purple")

png(filename ="figures/MvsF_DE_Top20_heatmap_decidua.png", width= 480, height= 580)
heatmap.2(logdge_CPM[i,], scale="row",
          labRow=y$genes$Geneid[i], labCol=y$targets$rep, 
          col=mycol, trace="none", density.info="none", 
          margin=c(6,6), lhei=c(2,10), dendrogram="column")
dev.off()
dev.off()

FvsM_topTable <- top.table$Geneid[1:100]
i <- which(y$genes$Geneid %in% FvsM_topTable)
mycol <- colorpanel(1000,"green","white","purple")

png(filename ="figures/MvsF_DE_Top100_heatmap_decidua.png", width= 480, height= 580)
heatmap.2(logdge_CPM[i,], scale="row",
          labRow=y$genes$Geneid[i], labCol=y$targets$rep, 
          col=mycol, trace="none", density.info="none", 
          margin=c(6,6), lhei=c(2,10), dendrogram="column")
dev.off()
dev.off()


#--------------------------------
# transcript violin jitter plots 
#--------------------------------
# Getting the sex of each sample from the metada
males <- pheno_ExRemovals[pheno_ExRemovals$offspring_sex == 'male',]
males <- unique(males["sample"])
males <- as.vector(unlist(males))
females <- pheno_ExRemovals[pheno_ExRemovals$offspring_sex == 'female',]
females <- unique(females["sample"])
females <- as.vector(unlist(females))

transformed <- y$E
transformedGenes<- y$genes
# Subsetting male and female data
transformed_males <- as.data.frame(transformed)[males]
transformed_females <- as.data.frame(transformed)[females]

# Adding columns with mean expression values
transformed_males$mean = apply(transformed_males,1,mean,na.rm=TRUE)
transformed_females$mean = apply(transformed_females,1,mean,na.rm=TRUE)

# New dataframe with gene names and male & female mean expression values
mean_exp <- data.frame(matrix(nrow=nrow(transformed)))
mean_exp$gene <- y$genes$gene
mean_exp$chr <- y$genes$chr
mean_exp$female <- transformed_females$mean
mean_exp$male <- transformed_males$mean
mean_exp_MvsF_All <- mean_exp 
t<-cbind(transformedGenes, transformed_females, transformed_males)

MAOA<- subset(t, gene_name=="MAOA")
TSIX<- subset(t, gene_name=="TSIX")

df_females = subset(MAOA, select = females)
df_males = subset(MAOA, select = males)
final_df_f <- as.data.frame(t(df_females))
final_df_m <- as.data.frame(t(df_males))
idk_m <- data.frame(final_df_m,sex=rep("male",ea=NROW(final_df_m)))
idk_f <- data.frame(final_df_f,sex=rep("female",ea=NROW(final_df_f)))

final_df<-rbind(idk_m,idk_f)
final_df$sex <- as.factor(final_df$sex)

colnames(final_df)[colnames(final_df)=="X200725"] <- "ENST00000542639.5"
colnames(final_df)[colnames(final_df)=="X200726"] <- "ENST00000338702.3"
colnames(final_df)[colnames(final_df)=="X200727"] <- "ENST00000497485.1"

means <- aggregate(ENST00000542639.5 ~  sex, final_df, mean)
var.test(ENST00000542639.5 ~ sex, final_df, alternative = "two.sided")
#res.ftest$p.value
#median <- aggregate(ENST00000542639.5 ~  sex, final_df, median)
png(filename ="figures/Violin_ENST00000542639.5_MAOA_decidua.png",width=6,height=6,units="in",res=1200)
p <- ggplot(final_df, aes(x=sex, y=ENST00000542639.5, color=sex))+ geom_violin() +scale_color_manual(values=c("black", "red", "gray"))
p + geom_boxplot(width=0.1, outlier.shape = NA) + geom_text(data = means, aes(label = round(ENST00000542639.5, digits = 2), y = ENST00000542639.5, color="black")) + stat_compare_means(method = "t.test") + geom_jitter(aes(shape = factor(sex)), size = 3, position=position_jitter(0.2)) + theme(legend.position="none")  + scale_shape_manual(values=c(14, 10)) + theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white")
)
dev.off()
dev.off()

means <- aggregate(ENST00000338702.3 ~  sex, final_df, mean)
res.ftest<-var.test(ENST00000338702.3 ~ sex, final_df, alternative = "two.sided")
res.ftest$p.value
#median <- aggregate(ENST00000338702.3 ~  sex, final_df, median)
png(filename ="figures/Violin_ENST00000338702.3_MAOA_decidua.png",width=6,height=6,units="in",res=1200)
p <- ggplot(final_df, aes(x=sex, y=ENST00000338702.3, color=sex))+ geom_violin() +scale_color_manual(values=c("black", "red", "gray"))
p + geom_boxplot(width=0.1, outlier.shape = NA) + geom_text(data = means, aes(label = round(ENST00000338702.3, digits = 2), y = ENST00000338702.3, color="black")) + stat_compare_means(method = "t.test") + geom_jitter(aes(shape = factor(sex)), size = 3, position=position_jitter(0.2)) + theme(legend.position="none")  + scale_shape_manual(values=c(14, 10)) + theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white")
)
dev.off()
dev.off()

means <- aggregate(ENST00000497485.1 ~  sex, final_df, mean)
res.ftest<-var.test(ENST00000497485.1 ~ sex, final_df, alternative = "two.sided")
res.ftest$p.value
#median <- aggregate(ENST00000338702.3 ~  sex, final_df, median)
png(filename ="figures/Violin_ENST00000497485.1_MAOA_decidua.png",width=6,height=6,units="in",res=1200)
p <- ggplot(final_df, aes(x=sex, y=ENST00000497485.1, color=sex))+ geom_violin() +scale_color_manual(values=c("black", "red", "gray"))
p + geom_boxplot(width=0.1, outlier.shape = NA) + geom_text(data = means, aes(label = round(ENST00000497485.1, digits = 2), y = ENST00000497485.1, color="black")) + stat_compare_means(method = "t.test") + geom_jitter(aes(shape = factor(sex)), size = 3, position=position_jitter(0.2)) + theme(legend.position="none")  + scale_shape_manual(values=c(14, 10)) + theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white")
)
dev.off()
dev.off()


  
TSIX<- subset(t, gene_name=="TSIX")

df_females = subset(TSIX, select = females)
df_males = subset(TSIX, select = males)
final_df_f <- as.data.frame(t(df_females))
final_df_m <- as.data.frame(t(df_males))
idk_m <- data.frame(final_df_m,sex=rep("male",ea=NROW(final_df_m)))
idk_f <- data.frame(final_df_f,sex=rep("female",ea=NROW(final_df_f)))

final_df<-rbind(idk_m,idk_f)
final_df$sex <- as.factor(final_df$sex)

colnames(final_df)[colnames(final_df)=="X202472"] <- "ENST00000604411.1"

means <- aggregate(ENST00000604411.1 ~  sex, final_df, mean)
var.test(ENST00000604411.1 ~ sex, final_df, alternative = "two.sided")
#res.ftest$p.value
#median <- aggregate(ENST00000542639.5 ~  sex, final_df, median)
png(filename ="figures/Violin_ENST00000604411.1_TSIX_decidua.png",width=6,height=6,units="in",res=1200)
p <- ggplot(final_df, aes(x=sex, y=ENST00000604411.1, color=sex))+ geom_violin() +scale_color_manual(values=c("black", "red", "gray"))
p + geom_boxplot(width=0.1, outlier.shape = NA) + geom_text(data = means, aes(label = round(ENST00000604411.1, digits = 2), y = ENST00000604411.1, color="black")) + stat_compare_means() + geom_jitter(aes(shape = factor(sex)), size = 3, position=position_jitter(0.2)) + theme(legend.position="none")  + scale_shape_manual(values=c(14, 10)) + theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white")
)
dev.off()
dev.off()

#--------------------------
# ARSD

ARSD<- subset(t, gene_name=="ARSD")

df_females = subset(ARSD, select = females)
df_males = subset(ARSD, select = males)
final_df_f <- as.data.frame(t(df_females))
final_df_m <- as.data.frame(t(df_males))
idk_m <- data.frame(final_df_m,sex=rep("male",ea=NROW(final_df_m)))
idk_f <- data.frame(final_df_f,sex=rep("female",ea=NROW(final_df_f)))

final_df<-rbind(idk_m,idk_f)
final_df$sex <- as.factor(final_df$sex)

colnames(final_df)[colnames(final_df)=="X199413"] <- "ENST00000381154.5"
colnames(final_df)[colnames(final_df)=="X199415"] <- "ENST00000495294.1"
colnames(final_df)[colnames(final_df)=="X199416"] <- "ENST00000495294.10"


means <- aggregate(ENST00000381154.5 ~  sex, final_df, mean)
var.test(ENST00000381154.5 ~ sex, final_df, alternative = "two.sided")
png(filename ="figures/Violin_ENST00000381154.5_ARSD_decidua.png",width=6,height=6,units="in",res=1200)
p <- ggplot(final_df, aes(x=sex, y=ENST00000381154.5, color=sex))+ geom_violin() +scale_color_manual(values=c("black", "red", "gray"))
p + geom_boxplot(width=0.1, outlier.shape = NA) + geom_text(data = means, aes(label = round(ENST00000381154.5, digits = 2), y = ENST00000381154.5, color="black")) + stat_compare_means(method = "t.test") + geom_jitter(aes(shape = factor(sex)), size = 3, position=position_jitter(0.2)) + theme(legend.position="none")  + scale_shape_manual(values=c(14, 10)) + theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white")
)
dev.off()
dev.off()

means <- aggregate(ENST00000495294.1 ~  sex, final_df, mean)
var.test(ENST00000495294.1 ~ sex, final_df, alternative = "two.sided")
png(filename ="figures/Violin_ENST00000495294.1_ARSD_decidua.png",width=6,height=6,units="in",res=1200)
p <- ggplot(final_df, aes(x=sex, y=ENST00000495294.1, color=sex))+ geom_violin() +scale_color_manual(values=c("black", "red", "gray"))
p + geom_boxplot(width=0.1, outlier.shape = NA) + geom_text(data = means, aes(label = round(ENST00000495294.1, digits = 2), y = ENST00000495294.1, color="black")) + stat_compare_means(method = "t.test") + geom_jitter(aes(shape = factor(sex)), size = 3, position=position_jitter(0.2)) + theme(legend.position="none")  + scale_shape_manual(values=c(14, 10)) + theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white")
)
dev.off()
dev.off()

means <- aggregate(ENST00000495294.10 ~  sex, final_df, mean)
var.test(ENST00000495294.10 ~ sex, final_df, alternative = "two.sided")
png(filename ="figures/Violin_EENST00000495294.10_ARSD_decidua.png",width=6,height=6,units="in",res=1200)
p <- ggplot(final_df, aes(x=sex, y=ENST00000495294.10, color=sex))+ geom_violin() +scale_color_manual(values=c("black", "red", "gray"))
p + geom_boxplot(width=0.1, outlier.shape = NA) + geom_text(data = means, aes(label = round(ENST00000495294.10, digits = 2), y = ENST00000495294.10, color="black")) + stat_compare_means(method="t.test") + geom_jitter(aes(shape = factor(sex)), size = 3, position=position_jitter(0.2)) + theme(legend.position="none")  + scale_shape_manual(values=c(14, 10)) + theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white")
)
dev.off()
dev.off()

