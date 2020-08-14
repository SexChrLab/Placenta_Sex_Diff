setwd(
  "~/Dropbox (ASU)/Placenta/BATCH2_PLACENTA_DECIDUA_ANALYSIS/HISAT_FeatureCounts/batch1_and_batch2/"
)
library(reshape)
library(ggpubr)
library(ggpubr)
library(dplyr)
library(lattice)
library(car)
library(splitstackshape)
library(ggrepel)

# read in data
# df_pvals_decidua_ImmuneDBgenes.txt
# df_pvals_term_ImmuneDBgenes.txt
innateFirstTri <- read.delim("innateDB/df_pvals_lateFirst_ImmuneDBgenes.txt")
names(innateFirstTri)[names(innateFirstTri) == "innateDB_inter_pisarska"] <- "Geneid"

innateTerm <- read.delim("innateDB/df_pvals_term_ImmuneDBgenes.txt")
names(innateTerm)[names(innateTerm) == "innateDB_inter_placenta"] <- "Geneid"

DEGsFirst <- read.delim("DEGs/pisarska/DEGs_placentas_genesFvsM_fdr1_lfc0.txt")
DEGsTerm <- read.delim("DEGs/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/DEGs_placentas_genesFvsM_fdr1_lfc0.txt")

innateTerm <- merge(innateTerm, DEGsTerm, by = "Geneid")
innateFirstTri <- merge(innateFirstTri, DEGsFirst, by = "Geneid")

#----- first trimester placentas 
innateFirstTri$FDRp <- p.adjust(innateFirstTri$wil_pval, n = length(innateFirstTri$wil_pval), method = "fdr")
sigFDR <- subset(innateFirstTri, FDRp <= 0.05)
dim(sigFDR)

up_female <- subset(innateFirstTri, logFC > .5)
up_female <- cbind(up_female, rep("FC =< 1", nrow(up_female)))
colnames(up_female)[20] <- "Color"

up_male <- subset(innateFirstTri, logFC < -.5)
up_male <- cbind(up_male, rep("FC >= 1", nrow(up_male)))
colnames(up_male)[20] <- "Color"

nonsig <- subset(innateFirstTri, logFC <.5 & logFC > -.5)
nonsig <- cbind(nonsig, rep("FC = 0", nrow(nonsig)))
colnames(nonsig)[20] <- "Color"

df_order <- rbind(up_female, up_male, nonsig)
colnames(df_order)[1] <- "Geneid"
df_order$Geneid <- factor(df_order$Geneid, levels = df_order$Geneid[order(df_order$logFC)])
df_order$Color <- as.factor(df_order$Color)
# innateFirstTri$color <- 

subset_data <- subset(df_order, logFC < -.5 | logFC > .5) 
subset_data_up <- subset(df_order, logFC > .5) 
subset_data_down <- subset(df_order, logFC < -.5) 

plot_first <- ggplot(df_order, aes(x=Geneid, y=logFC, color=Color)) + geom_point(size = 4) +
  scale_color_manual(values = c("#66C2A5", "#FC8D62","gray47")) + #"#66C2A5" greenish,  "#FC8D62" orange
  labs(title="first trimester female/male expression", x="gene ID", y="Fold change")+
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_blank()) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12)) +
  theme(axis.title=element_text(size=15)) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=15)) +
  theme(legend.position = "none") + scale_y_continuous(name="Fold change", limits=c(-3, 3), breaks=c(-3,-2.5, -2, -1.5, -1,-.5, 0, .5, 1, 1.5, 2, 2.5, 3)) + 
  geom_hline(yintercept = 0, colour="gray38", linetype="dashed", size = 1) + 
  geom_hline(yintercept = -.5, colour="#FC8D62", linetype="dashed", size = 1) + 
  geom_hline(yintercept = .5, colour="#66C2A5", linetype="dashed", size = 1) + 
  geom_text_repel(data=subset_data_up, aes(x = Geneid, y = logFC, label=Geneid), fontface="italic", segment.alpha = 1, box.padding = .5, color = "black", show.legend=FALSE, size = 4) +
  geom_text_repel(data=subset_data_down, aes(x = Geneid, y = logFC, label=Geneid), fontface="italic", segment.alpha = 1, box.padding = .5, color = "black", show.legend=FALSE, size = 4) +
  scale_x_discrete(expand = c(0.02,0.02)) +
annotate(geom="text", x=300, y=1.25, label="up-regulated 
female/male",
                           color="#66C2A5", size = 4.5) + 
annotate(geom="text", x=300, y=-1.5, label="down-regulated 
female/male",
                          color="#FC8D62", size = 4.5)
plot_first
ggsave("FC_placenta_firstTrimester_innateDB.pdf", plot_first,
       width = 5.5, height = 5.5)
dev.off()
dev.off()

#----- Term placentas 
innateTerm$FDRp <- p.adjust(innateTerm$wil_pval, n = length(innateTerm$wil_pval), method = "fdr")
sigFDR <- subset(innateTerm, FDRp <= 0.05)
dim(sigFDR)
innateTerm$logFC <- NULL
innateTerm$logFC <- log2((innateTerm$female_mean/innateTerm$male_mean))
up_female <- subset(innateTerm, logFC > .5)
up_female <- cbind(up_female, rep("FC =< 1", nrow(up_female)))
colnames(up_female)[14] <- "Color"

up_male <- subset(innateTerm, logFC < -.5)
up_male <- cbind(up_male, rep("FC >= 1", nrow(up_male)))
colnames(up_male)[14] <- "Color"

nonsig <- subset(innateTerm, logFC <.5 & logFC > -.5)
nonsig <- cbind(nonsig, rep("FC = 0", nrow(nonsig)))
colnames(nonsig)[14] <- "Color"

df_order <- rbind(up_female, up_male, nonsig)
colnames(df_order)[1] <- "Geneid"
df_order$Geneid <- factor(df_order$Geneid, levels = df_order$Geneid[order(df_order$logFC)])
df_order$Color <- as.factor(df_order$Color)
# innateTerm$color <- 

subset_data <- subset(df_order, logFC < -.5 | logFC > .5) 
subset_data_up <- subset(df_order, logFC > .5) 
subset_data_down <- subset(df_order, logFC < -.5) 

plot_term <- ggplot(df_order, aes(x=Geneid, y=logFC, color=Color)) + geom_point(size = 4) +
  scale_color_manual(values = c("#66C2A5", "#FC8D62", "gray47")) + #c("lightblue2", "lightpink2"))
  labs(title="term female/male expression", x="gene ID", y="Fold change")+
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_blank()) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12)) +
  theme(axis.title=element_text(size=15)) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=15)) +
  theme(legend.position = "none") + scale_y_continuous(name="Fold change", limits=c(-3, 3), breaks=c(-3, -2.5, -2, -1.5, -1,-.5, 0, .5, 1, 1.5, 2, 2.5, 3)) + 
  geom_hline(yintercept = 0, colour="gray38", linetype="dashed", size = 1) + 
  geom_hline(yintercept = -.5, colour="#FC8D62", linetype="dashed", size = 1) + 
  geom_hline(yintercept = .5, colour="#66C2A5", linetype="dashed", size = 1) + 
  geom_text_repel(data=subset_data_up, aes(x = Geneid, y = logFC, label=Geneid), fontface="italic", segment.alpha = 1, box.padding = .5, color = "black", show.legend=FALSE, size = 4) +
  geom_text_repel(data=subset_data_down, aes(x = Geneid, y = logFC, label=Geneid), fontface="italic", segment.alpha = 1, box.padding = .5, color = "black", show.legend=FALSE, size = 4) +
  scale_x_discrete(expand = c(0.02,0.02)) +
  annotate(geom="text", x=300, y=1.25, label="up-regulated 
female/male",
           color="#66C2A5", size = 4.5) + 
  annotate(geom="text", x=300, y=-1.5, label="down-regulated 
female/male",
           color="#FC8D62", size = 4.5)
plot_term
ggsave("FC_placenta_Term_innateDB.pdf", plot_term,
       width = 5.5, height = 5.5)
dev.off()
dev.off()

g <- ggarrange(plot_first, plot_term, ncol = 1, nrow = 2, align="h", 
          labels = c("A", "B"),
          common.legend = F)
ggsave("log2FC_innateDB_placentas.pdf", g)
dev.off()


