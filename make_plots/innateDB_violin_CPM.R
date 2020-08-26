setwd(
  "~/Dropbox (ASU)/Placenta/BATCH2_PLACENTA_DECIDUA_ANALYSIS/HISAT_FeatureCounts/batch1_and_batch2/"
)
library(reshape)
library(ggpubr)
library("gridExtra")

df_merged_pisarska <-
  read.delim ("genelists/pisarska/cpm_placenta_df_genes.txt")
df_merged_pisarska$group <- "late first trimester"
df_merged_placenta <-
  read.delim ("genelists/placentas_batch_1and2/cpm_placenta_df_genes.txt")
df_merged_placenta$group <- "term >= 36 weeks"
innateDB <- read.delim("innateDB/innateDB.txt", header = TRUE)

rbind.all.columns <- function(x, y) {
  x.diff <- setdiff(colnames(x), colnames(y))
  y.diff <- setdiff(colnames(y), colnames(x))
  x[, c(as.character(y.diff))] <- NA
  y[, c(as.character(x.diff))] <- NA
  return(rbind(x, y))
}

df_merged <- rbind.all.columns(df_merged_placenta, df_merged_pisarska)

innateDB_inter <- intersect(innateDB$Geneid, df_merged$Geneid)
genes <- c("XIST", "DDX3X", "CXCL10", "CXCL9", "EGR1", "MSR1", "JUN", "NLRC5", "PROCR", 
           "IL2RB", "SERPING1", "IL1R2", "TREML2", "CXCR4", "LEP", "AQP3", "PRKX", "MMP12")
#df_merged$value <- log10(df_merged$value)

#genes2 <- as.character(up_female$Geneid)
placenta_gene_PLOTS <- "./FIGURES/ViolinJitters/innateDB_genes_CPM.pdf"
pdf(placenta_gene_PLOTS)

# Function to plot histograms on top of each other
violoin_Func <- function(a) {
  geneDF <- subset(df_merged, Geneid == a)
  # means <- aggregate(value ~ sex, geneDF, mean)
  means <- aggregate(value ~group+sex, geneDF, mean)
  p <- ggplot(geneDF, aes(x = sex, y = value, color = sex)) +
    geom_violin() + scale_color_manual(values = c("black", "#66C2A5",  "#FC8D62")) +
    facet_wrap(~group) + theme(strip.text.x = element_text(size = 12)) +
    geom_boxplot(width = 0.1, outlier.shape = NA) + 
    geom_jitter(aes(shape = factor(sex)),
                size = 3,
                position = position_jitter(0.1)) +
    theme(legend.position = "none") + ggtitle(paste0(a, " cpm")) +
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
        label = round(value, digits = 2),
        y = value,
        color = "black"
      ),
      position = position_dodge(width = 0.9),
      size = 5
    ) +
    # theme(axis.text = element_text(size = 5, colour="black")) +
    stat_compare_means(method = "t.test",
                       label.x = 1.2,
                       label.y.npc = 1) +
    labs(y = "cpm")
}

# Map: iterates through items in Meta (which is a list of dataframes)
# and iterates through the names of the items in Meta simultaneously
violinPlots <- Map(violoin_Func, a = genes)
violinPlots
violinPlots$CXCL9
dev.off()
violinPlots$DDX3X


# DDX3X <- violinPlots$DDX3X
# JUN <- violinPlots$JUN
# MSR1 <- violinPlots$MSR1
# 
# grid.arrange(plot_first, DDX3X, EGR1, ncol = 2, nrow = 2, layout_matrix = rbind(c(1,1), c(2,3)))
