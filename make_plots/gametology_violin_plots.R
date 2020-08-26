setwd(
  "~/Dropbox (ASU)/Placenta/BATCH2_PLACENTA_DECIDUA_ANALYSIS/HISAT_FeatureCounts/batch1_and_batch2/"
)
library(reshape)
library(ggpubr)
library("gridExtra")

# read in placenta CPM data
df_merged_placenta <-
  read.delim ("genelists/placentas_batch_1and2/cpm_placenta_df_genes.txt")
df_merged_placenta$group <- "term >= 36 weeks"

# read in X and Y gametology gene ids
gametology <- read.delim("gametology/gametology.txt", header = TRUE)

# make an X and Y chromosome list of the gametology genes
gametology_inter_X <- intersect(gametology$X, df_merged_placenta$Geneid)
gametology_inter_Y <- intersect(gametology$Y, df_merged_placenta$Geneid)
# what is shared between the two lists 
gametology_inter <- unique(c(gametology_inter_X, gametology_inter_Y))

# subset the placent CPM data to only incldue the X and Y gametology genes
placenta_gametologys <- subset(df_merged_placenta, Geneid %in% gametology_inter)

#----------------
# for loop to format the placenta CPM into 
# the format needed for making plots
DF <- data.frame()
geneDF <- data.frame()
geneComb <- NULL
groupComb <- NULL 
variabledf <- NULL
Geneiddf <- NULL
sexdf <- NULL
CPM <- NULL


for(i in gametology_inter) {
  geneDF <- subset(placenta_gametologys, Geneid == i)
  variabledf <- c(variabledf, as.character(geneDF$variable))
  Geneiddf <- c(Geneiddf, as.character(geneDF$Geneid))
  sexdf <- c(sexdf, as.character(geneDF$sex))
  CPM <- c(CPM, as.numeric(geneDF$value))
  Y <- sub("X$", "Y", as.character(geneDF$Geneid))
  groupComb <- paste0(as.character(geneDF$Geneid),":",Y)
  geneComb <- c(geneComb, as.character(groupComb))
  DF <- cbind(variabledf, Geneiddf, sexdf, CPM, geneComb)
}

# save as a data frame
# rename the gene groups 
placenta_gametologys_df <- as.data.frame(DF)
names(placenta_gametologys_df)[names(placenta_gametologys_df) == "sexdf"] <- "sex"
names(placenta_gametologys_df)[names(placenta_gametologys_df) == "variabledf"] <- "sample"
names(placenta_gametologys_df)[names(placenta_gametologys_df) == "Geneiddf"] <- "gene"
placenta_gametologys_df$geneComb <- as.character(placenta_gametologys_df$geneComb)
placenta_gametologys_df[placenta_gametologys_df=="DDX3Y:DDX3Y"] <- "DDX3X:DDX3Y"
placenta_gametologys_df[placenta_gametologys_df=="PCDH11Y:PCDH11Y"] <- "PCDH11X:PCDH11Y"
placenta_gametologys_df[placenta_gametologys_df=="USP9Y:USP9Y"] <- "USP9X:USP9Y"
placenta_gametologys_df[placenta_gametologys_df=="ZFY:ZFY"] <- "ZFX:ZFY"
placenta_gametologys_df[placenta_gametologys_df=="KDM6A:UTY"] <- "UTX:UTY"
placenta_gametologys_df[placenta_gametologys_df=="KDM6A:KDM6A"] <- "UTX:UTY"
placenta_gametologys_df[placenta_gametologys_df=="UTY:UTY"] <- "UTX:UTY"

#----------------
# subset the placenta_gametologys_df by male and female 
male <- subset(placenta_gametologys_df, sex == "male")
female <- subset(placenta_gametologys_df, sex == "female")

# only X chromosome linked genes 
female_X <- subset(female, gene %in% gametology_inter_X)
male_X <- subset(male, gene %in% gametology_inter_X)

# only Y chromosome linked genes 
male_Y <- subset(male, gene %in% gametology_inter_Y)
female_Y <- subset(female, gene %in% gametology_inter_Y)

# Merge the male_X with male_Y
# this will be used to then sum the X and Y expression for each sample and each gene
maleXandY <- merge(x = male_X, y = male_Y, by = c("sample", "geneComb"), all = TRUE) #, by.y = names(geneComb))

# the X and Y gametology genes don't match since there is SRY and XIST
# this creates NAs
# replace NAs with zero 
maleXandY[is.na(maleXandY)] <- 0
# warning messages are okay

# sum the X and Y CPM values
maleXandY$CPM <- as.numeric(as.character(maleXandY$CPM.x)) + as.numeric(as.character(maleXandY$CPM.y))
# Drop the .y columns and the CPM.x 
drops <- c("CPM.x", "CPM.y","gene.y", "sex.y", "CPM.y")
maleDrops <- maleXandY[ , !(names(maleXandY) %in% drops)]

maledf <- as.data.frame(maleDrops)
# rename so of the .x columns 
male_rename <- rename(maledf, c("gene.x"="gene", "sex.x"="sex"))
# fill in missing male information 
male_rename$sex[is.na(male_rename$sex)] <- "male"
male_XY <- male_rename[complete.cases(male_rename[ , 3:4]),]

# add back in Y chromsome names
# Y chromosome gene names 
#male_rename_sry <- subset(male_rename, geneComb == "SRY:SRY")
#male_rename_sry$gene[is.na(male_rename_sry$gene)] <- "SRY"

# rbind male and female X chromosome expression information
male_female_Xchr <- rbind(male_XY, female_X)
# CPM column should be numeric
male_female_Xchr$CPM <- as.numeric(as.character(male_female_Xchr$CPM))
#--------------------------------------------
male_female_Xchr_PLOT <- "./FIGURES/male_female_Xchr.pdf"
pdf(male_female_Xchr_PLOT)

# Function to plot histograms on top of each other
violoin_Func <- function(a) {
  geneDF <- subset(male_female_Xchr, gene == a)
  means <- aggregate(CPM ~sex, geneDF, mean)
  p <- ggplot(geneDF, aes(x = sex, y = CPM, color = sex)) +
    geom_violin() + scale_color_manual(values = c("black", "#66C2A5",  "#FC8D62")) +
 #   facet_wrap(~geneComb) + theme(strip.text.x = element_text(size = 12)) +
    geom_boxplot(width = 0.1, outlier.shape = NA) + 
    geom_jitter(aes(shape = factor(sex)),
                size = 3,
                position = position_jitter(0.1)) +
    theme(legend.position = "none") + ggtitle(paste0(a, " CPM")) +
    theme(axis.title.x=element_text(size=15), 
          axis.text.x=element_text(size=12)) +
    theme(axis.title.y=element_text(size=15),
          axis.text.y=element_text(size=12)) +
    theme(axis.title=element_text(size=15)) +
    theme(legend.text=element_text(size=12)) +
    theme(legend.title=element_text(size=15)) +
    geom_text(
      data = means,
      aes(
        label = round(CPM, digits = 2),
        y = CPM,
        color = "black"
      ),
      position = position_dodge(width = 0.9),
      size = 5
    ) +
    theme(axis.text = element_text(size = 5, colour="black")) +
    stat_compare_means(method = "t.test",
                       label.x = 1.2,
                       label.y.npc = 1) +
    labs(y = "CPM")
}

# Map: iterates through items in Meta (which is a list of dataframes)
# and iterates through the names of the items in Meta simultaneously
violinPlots <- Map(violoin_Func, a = gametology_inter_X)
dev.off()
violinPlots$DDX3X
DDX3X <- violinPlots$DDX3X
PCDH11X <- violinPlots$PCDH11X
ZFX <- violinPlots$ZFX
USP9X <- violinPlots$USP9X
UTX <- violinPlots$KDM6A

ggsave("CPM_X+Y.pdf", grid.arrange(DDX3X, PCDH11X, USP9X, ZFX, UTX, ncol=5))
#---

# rbind male and female X chromosome expression information
male_female_Xchr_Only <- rbind(male_X, female_X)
# CPM column should be numeric
male_female_Xchr_Only$CPM <- as.numeric(as.character(male_female_Xchr_Only$CPM))
male_female_XchrOnly_PLOT <- "./FIGURES/male_female_Xchr.pdf"
pdf(male_female_XchrOnly_PLOT)

# Function to plot histograms on top of each other
violoin_Func <- function(a) {
  geneDF <- subset(male_female_Xchr_Only, gene == a)
  means <- aggregate(CPM ~sex, geneDF, mean)
  p <- ggplot(geneDF, aes(x = sex, y = CPM, color = sex)) +
    geom_violin() + scale_color_manual(values = c("black", "#66C2A5",  "#FC8D62")) +
    #   facet_wrap(~geneComb) + theme(strip.text.x = element_text(size = 12)) +
    geom_boxplot(width = 0.1, outlier.shape = NA) + 
    geom_jitter(aes(shape = factor(sex)),
                size = 3,
                position = position_jitter(0.1)) +
    theme(legend.position = "none") + ggtitle(paste0(a, " CPM")) +
    theme(axis.title.x=element_text(size=15), 
          axis.text.x=element_text(size=12)) +
    theme(axis.title.y=element_text(size=15),
          axis.text.y=element_text(size=12)) +
    theme(axis.title=element_text(size=15)) +
    theme(legend.text=element_text(size=12)) +
    theme(legend.title=element_text(size=15)) +
    geom_text(
      data = means,
      aes(
        label = round(CPM, digits = 2),
        y = CPM,
        color = "black"
      ),
      position = position_dodge(width = 0.9),
      size = 5
    ) +
    theme(axis.text = element_text(size = 5, colour="black")) +
    stat_compare_means(method = "t.test",
                       label.x = 1.2,
                       label.y.npc = 1) +
    labs(y = "CPM")
}

# Map: iterates through items in Meta (which is a list of dataframes)
# and iterates through the names of the items in Meta simultaneously
violinPlots <- Map(violoin_Func, a = gametology_inter_X)
dev.off()
violinPlots$DDX3X
DDX3X_x <- violinPlots$DDX3X
PCDH11X_x <- violinPlots$PCDH11X
ZFX_x <- violinPlots$ZFX
USP9X_x <- violinPlots$USP9X
UTX_x <- violinPlots$KDM6A


#grid.arrange(DDX3X, PCDH11X, USP9X, ZFX, UTX, ncol = 5, nrow = 5, layout_matrix = rbind(c(1,1), c(2,3)))
ggsave("CPM_XchrOnly.pdf", grid.arrange(DDX3X_x, PCDH11X_x, USP9X_x, ZFX_x, UTX_x, ncol=5, nrow = 2))
ggsave("CPM_X+Y.pdf", grid.arrange(DDX3X, PCDH11X, USP9X, ZFX, UTX, ncol=5))

ggsave("CPM_TopBottom.pdf", grid.arrange(DDX3X_x, PCDH11X_x, USP9X_x, ZFX_x, UTX_x, DDX3X, PCDH11X, USP9X, ZFX, UTX, ncol=5, nrow = 2))

#---------------------------------------------------------
