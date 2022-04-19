setwd("~/Dropbox (ASU)/Placenta/BATCH2_PLACENTA_DECIDUA_ANALYSIS/HISAT_FeatureCounts/batch1_and_batch2/")
library(VennDiagram)
library(RColorBrewer)

PisarskaALL <- read.table("DEGs/pisarska/DEGs_placentas_genes_fdr1_lfc0_NOFPKMFILTER.txt", header = TRUE)
WilsonALL <- read.table("DEGs/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/DEGs_placentas_genes_fdr1_lfc0_NOFPKMFILTER.txt", header = TRUE)

PisarskaDEG <- read.table("DEGs/pisarska/DEGs_placentas_genesFvsM_fdr05_lfc0.txt", header = TRUE)
PisarskaDEG_genename <- PisarskaDEG$Geneid
WilsonDEG <- read.table("DEGs/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/DEGs_placentas_genesFvsM_fdr05_lfc0.txt", header = TRUE)
WilsonDEG_genename <- WilsonDEG$Geneid
uniquePisarskaDEG_genename <- setdiff(PisarskaDEG_genename, WilsonDEG_genename)
uniqueWilsonDEG_genename <- setdiff(WilsonDEG_genename, PisarskaDEG_genename)
sharedGenes <- intersect(WilsonDEG_genename, PisarskaDEG_genename)

df <- merge(WilsonALL, PisarskaALL, by = "Geneid")
df.Y <- subset(df, chr.x == "chrY") #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(1, nrow(df.Y)))
colnames(df.Y)[18] <- "color"
df.X <- subset(df, chr.x == "chrX") #define Y-chromosomal, blue
df.X <- cbind(df.X, rep(2, nrow(df.X)))
colnames(df.X)[18] <- "color"
df.A <- subset(df, chr.x != "chrX" & chr.x != "chrY" ) #define Y-chromosomal, blue
df.A <- cbind(df.A, rep(3, nrow(df.A)))
colnames(df.A)[18] <- "color"
df.t <- rbind(df.Y, df.X, df.A)
df.t$color <- as.factor(df.t$color)
             
shared <- df.t[df.t$Geneid %in% sharedGenes, ]
uniqueToWilson <- df.t[df.t$Geneid %in% uniqueWilsonDEG_genename, ]
uniqueToGonzalez <- df.t[df.t$Geneid %in% uniquePisarskaDEG_genename, ]

png(
  filename = "Venn_Correlation_DEGs_datasets/log2FC_correlation_Gonzalez_vs_Wilson_sharedFDR05_r.png",
  width = 4,
  height = 4,
  units = "in",
  res = 1200
)

p <- ggplot(data = shared, aes(x = logFC.x, y = logFC.y, color=color)) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "gray29",
    formula = y ~ x,
    size = .5
  ) +
  geom_point(size = 4) + xlim(-13.5, 8.5) + ylim(-13.5, 8.5) +   
  scale_color_manual(values = c("#D95F02", "#66A61E", "#7570B3"))
p
lm_eqn <- function(shared) {
  m <- lm(logFC.y ~ logFC.x, shared)
  eq <- list(r2 = format(summary(m)$r.squared, digits = 3))
  as.character(as.expression(eq))
}
lm_eqn(shared)
r2 <- lm_eqn(shared)
p1 <-
  p + annotate(
    "text",
    x = -6,
    y = 3,
    label = paste("~italic(r)==", r2),
    parse = TRUE,
    color = "gray29",
    size = 7
  )
p1 + labs(
  title = "",
  x = expression(log[2](FC) ~ "term"),
  y = expression(log[2](FC) ~ "first trimester")
) +
  theme(plot.title = element_text(size = 18)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 18))
dev.off()
dev.off()

png(
  filename = "Venn_Correlation_DEGs_datasets/log2FC_correlation_Gonzalez_vs_Wilson_uniqueToGonzalezFDR05_r.png",
  width = 4,
  height = 4,
  units = "in",
  res = 1200
)

p <- ggplot(data = uniqueToGonzalez, aes(x = logFC.x, y = logFC.y, color = color)) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "gray29",
    formula = y ~ x,
    size = .5
  ) +
  geom_point(size = 4) + xlim(-13.5, 8.5) + ylim(-13.5, 8.5) +
  scale_color_manual(values = c("#D95F02", "#66A61E", "#7570B3"))
lm_eqn <- function(uniqueToGonzalez) {
  m <- lm(logFC.y ~ logFC.x, uniqueToGonzalez)
  eq <- list(r2 = format(summary(m)$r.squared, digits = 3))
  as.character(as.expression(eq))
}
lm_eqn(uniqueToGonzalez)
r2 <- lm_eqn(uniqueToGonzalez)
p1 <-
  p + annotate(
    "text",
    x = -6,
    y = 3,
    label = paste("~italic(r)==", r2),
    parse = TRUE,
    color = "gray29",
    size = 7
  )
p1 + labs(
  title = "",
  x = expression(log[2](FC) ~ "term"),
  y = expression(log[2](FC) ~ "first trimester")
) +
  theme(plot.title = element_text(size = 18)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 18))
dev.off()
dev.off()

png(
  filename = "Venn_Correlation_DEGs_datasets/log2FC_correlation_Gonzalez_vs_Wilson_uniqueToWilsonFDR05_r.png",
  width = 4,
  height = 4,
  units = "in",
  res = 1200
)

p <- ggplot(data = uniqueToWilson, aes(x = logFC.x, y = logFC.y, color = color)) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "gray29",
    formula = y ~ x,
    size = .5
  ) +
  geom_point(size = 4) + xlim(-13.5, 8.5) + ylim(-13.5, 8.5) +
  scale_color_manual(values = c("#D95F02", "#66A61E", "#7570B3"))
lm_eqn <- function(uniqueToWilson) {
  m <- lm(logFC.y ~ logFC.x, uniqueToWilson)
  eq <- list(r2 = format(summary(m)$r.squared, digits = 3))
  as.character(as.expression(eq))
}
lm_eqn(uniqueToWilson)
r2 <- lm_eqn(uniqueToWilson)
p1 <-
  p + annotate(
    "text",
    x = -6,
    y = 3,
    label = paste("~italic(r)==", r2),
    parse = TRUE,
    color = "gray29",
    size = 7
  )
p1 + labs(
  title = "",
  x = expression(log[2](FC) ~ "term"),
  y = expression(log[2](FC) ~ "first trimester")
) +
  theme(plot.title = element_text(size = 18)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 18))
dev.off()
dev.off()

