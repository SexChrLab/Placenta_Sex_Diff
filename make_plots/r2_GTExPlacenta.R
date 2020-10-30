# Sex differences in placenta compared to sex differences in adult GTEx tissues.
# Lopes-Ramos et al. 2020 gene lists per tissue

# Set working directory
setwd("~/Dropbox (ASU)/Placenta/BATCH2_PLACENTA_DECIDUA_ANALYSIS/HISAT_FeatureCounts/batch1_and_batch2/sexDiff_GTEx_version6/")

# load libraries
library(dplyr)
library(tidyverse)
#install.packages("corrplot")
source("http://www.sthda.com/upload/rquery_cormat.r")
library(corrplot)
library(UpSetR)


# Notes from Lopes-Ramos et al. 2020 paper
# 29 tissues analyzes. Version 6 of GTEx
# Expression information for 30,243 genes, including 1,074 encoded on the sex chromosomes
# Only considered tissues with samples collected from both men and women
# More than 30 samples in either of the sexes. This included 8,279 samples representing 29 tissue types from 548 research subjects (360 males and 188 females).

# voom for differential expression
# efined differential expression using a cutoff
# of FDR < 0.05 and absolute fold change >= 1.5

# Even though most transcription factors (TFs) are not differentially expressed between males and females, many have sex-biased regulatory targeting patterns
# Compared the networks between males and females in each tissue and found significant sex differences in gene regu- latory networks across all tissues

## To study sex bias in regulatory processes
# 1) performed differential expression analysis
# 2) followed by gene regulatory network analysis

#------
# 1) sex-bias expression across tissues
#------
# read in the each GTEx tissue file
# each file will have the geneID "ID" and the log2M:F expression "logFC"
GTEx_Tissues <-
  c(
    "Whole_blood",
    "Thyroid",
    "Stomach",
    "Spleen",
    "Intestine_terminal_ileum",
    "Skin",
    "Pituitary",
    "Pancreas",
    "Tibial_nerve",
    "Skeletal_muscle",
    "Lung",
    "Liver",
    "Heart_left_ventricle",
    "Heart_atrial_appendage",
    "Esophagus_muscularis",
    "Esophagus_musosa",
    "Gastroesophageal_junction",
    "Colon_tranverse",
    "Colon_sigmoid",
    "Breast",
    "Brain_basal_ganglia",
    "Brain_cerebellum",
    "Brain_other",
    "Artery_tibial",
    "Artery_coronary",
    "Artery_aorta",
    "Adrenal_gland",
    "Adipose_visceral",
    "Adipose_subcutaneous"
  )


GTEx_Tissues <-
  c(
    "Whole_blood",
    "Brain_basal_ganglia",
    "Brain_cerebellum",
    "Brain_other"
  )

for (i in GTEx_Tissues) {
  tissue <- read.delim(paste0(i, ".txt"))
  tissue <- tissue[2:3] # select columns 2 and 3, the ID and logFC
  assign(paste(i, sep = ''), tissue) # assign each tissue as a dataframe in r
  # variable name for each tissue data frame is defined in i, GTEx_Tissues
}
# read in the placenta data
Placenta <- read.delim("Placenta.txt")
Placenta_05 <- subset(Placenta, Placenta$adj.p.value <= 0.05)

# only select the columns of interest, ID and logFC
Placenta <- Placenta[1:2]
# the placenta data is log2M:F and the GTEx is log2F:M
# multiple the logFC by -1  to be in the same direction as the GTEx data
Placenta$logFC <- (Placenta$logFC * -1) 

# Merge all tissue files by geneID
# if a tissue file doesn't have that geneID (i.e that gene is NOT sex differentially expressed in that tissue)
# than for that gene put NA
MergeTissuesByID <-
  list(
    Placenta,
    Whole_blood,
    Thyroid,
    Stomach,
    Spleen,
    Intestine_terminal_ileum,
    Skin,
    Pituitary,
    Pancreas,
    Tibial_nerve,
    Skeletal_muscle,
    Lung,
    Liver,
    Heart_left_ventricle,
    Heart_atrial_appendage,
    Esophagus_muscularis,
    Esophagus_musosa,
    Gastroesophageal_junction,
    Colon_tranverse,
    Colon_sigmoid,
    Breast,
    Brain_basal_ganglia,
    Brain_cerebellum,
    Brain_other,
    Artery_tibial,
    Artery_coronary,
    Artery_aorta,
    Adrenal_gland,
    Adipose_visceral,
    Adipose_subcutaneous
  ) %>% reduce(full_join, by = "ID")

MergeTissuesByID <-
  list(
    Placenta,
    Whole_blood,
    Brain_basal_ganglia,
    Brain_cerebellum,
    Brain_other
  ) %>% reduce(full_join, by = "ID")




# set the row names to be the gene IDs
rownames(MergeTissuesByID) <- MergeTissuesByID[, 1]

# remove all of the columsn that contain the geneID information
# the geneIDs are no longer needed since this information is saved
# as the row names
df_tissues <-
  MergeTissuesByID[,-grep("ID", colnames(MergeTissuesByID))]
# rename the columns to correspond to the tissue that column information belongs to
colnames(df_tissues) <-
  c(
    "Placenta",
    "Whole_blood",
    "Thyroid",
    "Stomach",
    "Spleen",
    "Intestine_terminal_ileum",
    "Skin",
    "Pituitary",
    "Pancreas",
    "Tibial_nerve",
    "Skeletal_muscle",
    "Lung",
    "Liver",
    "Heart_left_ventricle",
    "Heart_atrial_appendage",
    "Esophagus_muscularis",
    "Esophagus_musosa",
    "Gastroesophageal_junction",
    "Colon_tranverse",
    "Colon_sigmoid",
    "Breast",
    "Brain_basal_ganglia",
    "Brain_cerebellum",
    "Brain_other",
    "Artery_tibial",
    "Artery_coronary",
    "Artery_aorta",
    "Adrenal_gland",
    "Adipose_visceral",
    "Adipose_subcutaneous"
  )

colnames(df_tissues) <-
  c(
    "Placenta",
    "Whole_blood",
    "Brain_basal_ganglia",
    "Brain_cerebellum",
    "Brain_other"
  )



# correlation and plot
#rquery.cormat(df_tissues, type = "full")
rquery.cormat(df_tissues)



#-- correlation matrix with corrplot 
tissueMatrix <- as.matrix(df_tissues)
corrplot(tissueMatrix, is.corr = FALSE)
corrplot(cor(tissueMatrix)) #, is.corr = FALSE)

col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                           "cyan", "#007FFF", "blue", "#00007F"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
col3 <- colorRampPalette(c("red", "white", "blue")) 
col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F",
                           "cyan", "#007FFF", "blue", "#00007F"))
whiteblack <- c("white", "black")

## using these color spectra
corrplot(cor(na.omit(tissueMatrix)), order = "hclust", addrect = 5, col = col1(100))
corrplot(cor(na.omit(tissueMatrix)), order = "hclust", addrect = 2, col = col3(20))
corrplot(cor(na.omit(tissueMatrix)), order = "hclust", addrect = 2, col = heat.colors(100))
corrplot(cor(na.omit(tissueMatrix)), order = "hclust", addrect = 2,  col = col4(20), method = "circle", cl.lim = c(0.89, 1))

p.mat <- cor.mtest(cor(na.omit(tissueMatrix)))$p
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(cor(na.omit(tissueMatrix)), is.corr=FALSE, method = "color", col = col(200),
         cl.lim = c(0.90, 1),
         type = "upper", order = "hclust", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         number.digits = 4, 
        # addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
       #  p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag = FALSE)

#--- UpSet plot

listInput <- list(Placenta = Placenta_05$ID,
                  Whole_blood = Whole_blood$ID,
                  Thyroid = Thyroid$ID,
                  Stomach = Stomach$ID,
                  Spleen = Spleen$ID,
                  Intestine_terminal_ileum = Intestine_terminal_ileum$ID,
                  Skin = Skin$ID,
                  Pituitary = Pituitary$ID,
                  Pancreas = Pancreas$ID,
                  Tibial_nerve = Tibial_nerve$ID,
                  Skeletal_muscle = Skeletal_muscle$ID,
                  Lung = Lung$ID,
                  Liver = Liver$ID,
                  Heart_left_ventricle = Heart_left_ventricle$ID,
                  Heart_atrial_appendage = Heart_atrial_appendage$ID,
                  Esophagus_muscularis = Esophagus_muscularis$ID,
                  Esophagus_musosa = Esophagus_musosa$ID,
                  Gastroesophageal_junction = Gastroesophageal_junction$ID,
                  Colon_tranverse = Colon_tranverse$ID,
                  Colon_sigmoid = Colon_sigmoid$ID,
                  # Breast = Breast$ID,
                  Brain_basal_ganglia = Brain_basal_ganglia$ID,
                  Brain_cerebellum = Brain_cerebellum$ID,
                  Brain_other = Brain_other$ID,
                  Artery_tibial = Artery_tibial$ID,
                  Artery_coronary = Artery_coronary$ID,
                  Artery_aorta = Artery_aorta$ID,
                  Adrenal_gland = Adrenal_gland$ID,
                  Adipose_visceral = Adipose_visceral$ID,
                  Adipose_subcutaneous = Adipose_subcutaneous$ID)

upset(fromList(listInput), point.size = .5, line.size = 1, 
      order.by = "freq", nsets = 20, 
      number.angles = 0, text.scale = c(2, 2, 2, 1, 2, 1.5)) #order.by = "freq", 



listInput <- list(Placenta = Placenta$ID,
                  Whole_blood = Whole_blood$ID,
                  Brain_basal_ganglia = Brain_basal_ganglia$ID,
                  Brain_cerebellum = Brain_cerebellum$ID,
                  Brain_other = Brain_other$ID)
#pdf(file="Clark_DEG_vs_ASE_withinEachHybrid.pdf")
upset(fromList(listInput), point.size = 2.5, line.size = 1, 
      order.by = "freq", nsets = 10, 
      number.angles = 0, text.scale = c(2, 2, 2, 1, 2, 1.5)) #order.by = "freq", 
#dev.off()
#dev.off()
intersect(Placenta$ID, )

