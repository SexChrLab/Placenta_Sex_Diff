library(grid)
library(gridExtra)
library(ggplot2)
setwd("~/Dropbox (ASU)/Placenta/BATCH2_PLACENTA_DECIDUA_ANALYSIS/HISAT_FeatureCounts/batch1_and_batch2/")

# read in placenta DEGs files that we gerenated
Wilson <-
  read.delim(
    "DEGs/placentas_batch_1and2/batch_BirthWeight_lane_PC1_PC2/DEGs_placentas_genesFvsM_fdr1_lfc0.txt"
  )
Gonzalez <-
  read.delim("DEGs/pisarska/DEGs_placentas_genesFvsM_fdr1_lfc0.txt")

# read in the Gong et al. 2012 gene lists with the log2M:F information
# Supplemental table 2: Differentially expressed genes (DEGs) in male and female placentas
Gong <- read.delim("GongEtAl_2018/GongEtAl_SupTab2.txt")

# Supplemental table 3: X-chromosome genes with sex-biased expression in 19 GTEx tissues and placenta
#Sup3_Gong <- read.delim("GongEtAl_2012/GongEtAl_SupTab3.txt")

# Supplemental table 4: Expression levels of the 22 genes with female-biased expression uniquely in placenta
#Sup4_Gong <- read.delim("GongEtAl_2012/GongEtAl_SupTab4.txt")

# A) r2 for genes in Gong et al. Supplemental table 2: Differentially expressed genes (DEGs) in male and female placentas
# NOTE Gong et al used an adj p-value < 0.01 threshold.
# The DEGs for the first trimester and full-term described here
# do not have p-value or  log2FC cut off.

# If the gene in the Gong et al paper is also present in the first or full
# term placentas, then keep that gene
GongGonzalez <- intersect(Gong$Geneid, Gonzalez$Geneid)
GongWilson <- intersect(Gong$Geneid, Wilson$Geneid)

# Wilson and Gong
GongWilson <- subset(Wilson, Wilson$Geneid %in% GongWilson)
GongGonzalez <-
  subset(Gonzalez, Gonzalez$Geneid %in% GongGonzalez)

df <- merge(GongWilson, GongGonzalez, by = "Geneid")
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

p <- ggplot(data = df.t, aes(x = logFC.x, y = logFC.y, color=color)) + # Wilson is logFC.x 
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "gray29",
    formula = y ~ x,
    size = .5
  ) +
  geom_point(size = 4) + xlim(-2, 8.5) + ylim(-2, 8.5) +
  scale_color_manual(values = c("#66A61E","#7570B3"))#"#D95F02",
p
lm_eqn <- function(df.t) {
  m <- lm(logFC.y ~ logFC.x, df.t)
  eq <- list(r2 = format(summary(m)$r.squared, digits = 3))
  as.character(as.expression(eq))
}
lm_eqn(df.t)
r2 <- lm_eqn(df.t)
p1 <-
  p + annotate(
    "text",
    x = 2,
    y = 7,
    label = paste("~italic(r)^2==", r2),
    parse = TRUE,
    color = "gray29",
    size = 5
  )
plot1 <- p1 + labs(
  title = "Wilson vs. Gonzalez",
  x = expression(log[2](FC) ~ "Gonzalez"),
  y = expression(log[2](FC) ~ "Wilson")
) +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))
plot1
ggsave("FIGURES/Venn_Correlation_DEGs_datasets/log2FC_correlation_Wilson_vs_Gonsalex_Sup2_SDE.pdf", plot1, width = 3.5, height = 3.5)

df <- merge(GongWilson, Gong, by = "Geneid")
df.Y <- subset(df, chr == "chrY") #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(1, nrow(df.Y)))
colnames(df.Y)[11] <- "color"
df.X <- subset(df, chr == "chrX") #define Y-chromosomal, blue
df.X <- cbind(df.X, rep(2, nrow(df.X)))
colnames(df.X)[11] <- "color"
df.A <- subset(df, chr != "chrX" & chr != "chrY" ) #define Y-chromosomal, blue
df.A <- cbind(df.A, rep(3, nrow(df.A)))
colnames(df.A)[11] <- "color"
df.t <- rbind(df.Y, df.X, df.A)
df.t$color <- as.factor(df.t$color)

p <- ggplot(data = df.t, aes(x = logFC.x, y = logFC.y, color=color)) + # Gonzalez is logFC.x 
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "gray29",
    formula = y ~ x,
    size = .5
  ) +
  geom_point(size = 4) + xlim(-2, 8.5) + ylim(-2, 8.5) +
  scale_color_manual(values = c("#66A61E","#7570B3"))#"#D95F02",
p
lm_eqn <- function(df.t) {
  m <- lm(logFC.y ~ logFC.x, df.t)
  eq <- list(r2 = format(summary(m)$r.squared, digits = 3))
  as.character(as.expression(eq))
}
lm_eqn(df.t)
r2 <- lm_eqn(df.t)
p2 <-
  p + annotate(
    "text",
    x = 2,
    y = 7,
    label = paste("~italic(r)^2==", r2),
    parse = TRUE,
    color = "gray29",
    size = 5
  )
plot2 <- p2 + labs(
  title = "Wilson vs. Gong",
  x = expression(log[2](FC) ~ "Gong"),
  y = expression(log[2](FC) ~ "Wilson")
) +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))
plot2
ggsave("FIGURES/Venn_Correlation_DEGs_datasets/log2FC_correlation_Wilson_vs_Gong_Sup2_SDE.pdf", plot2, width = 3.5, height = 3.5)


#  Gonzalez vs Gong
df <- merge(GongGonzalez, Gong, by = "Geneid")
df.Y <- subset(df, chr == "chrY") #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(1, nrow(df.Y)))
colnames(df.Y)[11] <- "color"
df.X <- subset(df, chr == "chrX") #define Y-chromosomal, blue
df.X <- cbind(df.X, rep(2, nrow(df.X)))
colnames(df.X)[11] <- "color"
df.A <- subset(df, chr != "chrX" & chr != "chrY" ) #define Y-chromosomal, blue
df.A <- cbind(df.A, rep(3, nrow(df.A)))
colnames(df.A)[11] <- "color"
df.t <- rbind(df.Y, df.X, df.A)
df.t$color <- as.factor(df.t$color)

p <- ggplot(data = df.t, aes(x = logFC.x, y = logFC.y, color=color)) + # Wilson is logFC.x 
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "gray29",
    formula = y ~ x,
    size = .5
  ) +
  geom_point(size = 4) + xlim(-2, 8.5) + ylim(-2, 8.5) +
  scale_color_manual(values = c("#66A61E","#7570B3"))#"#D95F02",
p
lm_eqn <- function(df.t) {
  m <- lm(logFC.y ~ logFC.x, df.t)
  
  eq <- list(r2 = format(summary(m)$r.squared, digits = 3))
  as.character(as.expression(eq))
  
}
lm_eqn(df.t)
r2 <- lm_eqn(df.t)
p3 <-
  p + annotate(
    "text",
    x = 2,
    y = 7,
    label = paste("~italic(r)^2==", r2),
    parse = TRUE,
    color = "gray29",
    size = 5 
  )
plot3 <- p3 + labs(
  title = "Gonzalez vs. Gong",
  x = expression(log[2](FC) ~ "Gong"),
  y = expression(log[2](FC) ~ "Gonzalez")
) +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))
plot3
ggsave("FIGURES/Venn_Correlation_DEGs_datasets/log2FC_correlation_Gonzalez_vs_Gong_Sup3_SDE.pdf", plot3, width = 3.5, height = 3.5)

# ggsave("filename.pdf", 
#        grid.arrange(plot1, plot2, plot3, nrow=1, widths=c(3, 3)),
#        width = 15, height = 15)

# grid.arrange(
#   plot1, plot2, plot3,
#   widths = c(2, 2, 2),
#   layout_matrix = rbind(c(1, 2, 3))
# )
# ggsave("filename.pdf",
#        grid.arrange(
#          plot1, plot2, plot3,
#          widths = c(2, 2, 2),
#          layout_matrix = rbind(c(1, 2, 3)),
#          widht = 11, height = 3))

lay = rbind(c(1, 2, 3))
ggsave("filename.pdf", 
       grid.arrange(grobs = list(plot1, plot2, plot3), 
             layout_matrix = lay),
       width = 7.5, height = 2.67)

#---------------------------------------------------------------------------
#--- Gonzalez and Gong
Gonzalez_Sup2 <- subset(Gonzalez, Gonzalez$Geneid %in% Sup2_GongGonzalez)
Gong_Sup2 <-
  subset(Sup2_Gong, Sup2_Gong$Geneid %in% Sup2_GongGonzalez)

df <- merge(Gonzalez_Sup2, Gong_Sup2, by = "Geneid")
df.Y <- subset(df, chr == "chrY") #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(1, nrow(df.Y)))
colnames(df.Y)[11] <- "color"
df.X <- subset(df, chr == "chrX") #define Y-chromosomal, blue
df.X <- cbind(df.X, rep(2, nrow(df.X)))
colnames(df.X)[11] <- "color"
df.A <- subset(df, chr != "chrX" & chr != "chrY" ) #define Y-chromosomal, blue
df.A <- cbind(df.A, rep(3, nrow(df.A)))
colnames(df.A)[11] <- "color"
df.t <- rbind(df.Y, df.X, df.A)
df.t$color <- as.factor(df.t$color)


png(
  filename = "FIGURES/Venn_Correlation_DEGs_datasets/log2FC_correlation_Gong_vs_Gonzalez_Sup2_SDE.png",
  width = 4,
  height = 4,
  units = "in",
  res = 1200
)

p <- ggplot(data = df.t, aes(x = logFC.x, y = logFC.y, color=color)) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "gray29",
    formula = y ~ x,
    size = .5
  ) +
  geom_point(size = 4) + xlim(-2, 8.5) + ylim(-2, 8.5) +
  scale_color_manual(values = c("#66A61E","#7570B3"))#"#D95F02",
p
lm_eqn <- function(df.t) {
  m <- lm(logFC.y ~ logFC.x, df.t)
  
  eq <- list(r2 = format(summary(m)$r.squared, digits = 3))
  as.character(as.expression(eq))
  
}
lm_eqn(df.t)
r2 <- lm_eqn(df.t)
p1 <-
  p + annotate(
    "text",
    x = 3,
    y = 5,
    label = paste("~italic(r)^2==", r2),
    parse = TRUE,
    color = "gray29",
    size = 7
  )
p1 + labs(
  title = "Gong vs. Gonzalez",
  x = expression(log[2](FC) ~ "Gonzalez"),
  y = expression(log[2](FC) ~ "Gong")
) +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))

dev.off()
dev.off()

# If the gene in the Gong et al paper is also present in the first or full
# term placentas, then keep that gene
Sup3_GongGonzalez <- intersect(Sup3_Gong$Geneid, Gonzalez$Geneid)
Sup3_GongWilson <- intersect(Sup3_Gong$Geneid, Wilson$Geneid)

# Wilson and Gong
Wilson_Sup3 <- subset(Wilson, Wilson$Geneid %in% Sup3_GongWilson)
Gong_Sup3 <-
  subset(Sup3_Gong, Sup3_Gong$Geneid %in% Sup3_GongWilson)

df <- merge(Wilson_Sup3, Gong_Sup3, by = "Geneid")
df.Y <- subset(df, chr == "chrY") #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(1, nrow(df.Y)))
colnames(df.Y)[11] <- "color"
df.X <- subset(df, chr == "chrX") #define Y-chromosomal, blue
df.X <- cbind(df.X, rep(2, nrow(df.X)))
colnames(df.X)[11] <- "color"
df.A <- subset(df, chr != "chrX" & chr != "chrY" ) #define Y-chromosomal, blue
df.A <- cbind(df.A, rep(3, nrow(df.A)))
colnames(df.A)[11] <- "color"
df.t <- rbind(df.Y, df.X, df.A)
df.t$color <- as.factor(df.t$color)


png(
  filename = "FIGURES/Venn_Correlation_DEGs_datasets/log2FC_correlation_Gong_vs_Wilson_Sup3_SDE.png",
  width = 4,
  height = 4,
  units = "in",
  res = 1200
)

p <- ggplot(data = df.t, aes(x = logFC.x, y = logFC.y, color=color)) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "gray29",
    formula = y ~ x,
    size = .5
  ) +
  geom_point(size = 4) + xlim(-2, 8.5) + ylim(-2, 8.5) +
  scale_color_manual(values = c("#66A61E","#7570B3"))#"#D95F02",
p
lm_eqn <- function(df.t) {
  m <- lm(logFC.y ~ logFC.x, df.t)
  
  eq <- list(r2 = format(summary(m)$r.squared, digits = 3))
  as.character(as.expression(eq))
  
}
lm_eqn(df.t)
r2 <- lm_eqn(df.t)
p1 <-
  p + annotate(
    "text",
    x = 3,
    y = 5,
    label = paste("~italic(r)^2==", r2),
    parse = TRUE,
    color = "gray29",
    size = 7
  )
p1 + labs(
  title = "Gong vs. Wilson",
  x = expression(log[2](FC) ~ "Wilson"),
  y = expression(log[2](FC) ~ "Gong")
) +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))

dev.off()
dev.off()

#--- Gonzalez and Gong
Gonzalez_Sup3 <- subset(Gonzalez, Gonzalez$Geneid %in% Sup3_GongGonzalez)
Gong_Sup3 <-
  subset(Sup3_Gong, Sup3_Gong$Geneid %in% Sup3_GongGonzalez)

df <- merge(Gonzalez_Sup3, Gong_Sup3, by = "Geneid")
df.Y <- subset(df, chr == "chrY") #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(1, nrow(df.Y)))
colnames(df.Y)[11] <- "color"
df.X <- subset(df, chr == "chrX") #define Y-chromosomal, blue
df.X <- cbind(df.X, rep(2, nrow(df.X)))
colnames(df.X)[11] <- "color"
df.A <- subset(df, chr != "chrX" & chr != "chrY" ) #define Y-chromosomal, blue
df.A <- cbind(df.A, rep(3, nrow(df.A)))
colnames(df.A)[11] <- "color"
df.t <- rbind(df.Y, df.X, df.A)
df.t$color <- as.factor(df.t$color)


png(
  filename = "FIGURES/Venn_Correlation_DEGs_datasets/log2FC_correlation_Gong_vs_Gonzalez_Sup3_SDE.png",
  width = 4,
  height = 4,
  units = "in",
  res = 1200
)

p <- ggplot(data = df.t, aes(x = logFC.x, y = logFC.y, color=color)) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "gray29",
    formula = y ~ x,
    size = .5
  ) +
  geom_point(size = 4) + xlim(-2, 8.5) + ylim(-2, 8.5) +
  scale_color_manual(values = c("#66A61E","#7570B3"))#"#D95F02",
p
lm_eqn <- function(df.t) {
  m <- lm(logFC.y ~ logFC.x, df.t)
  
  eq <- list(r2 = format(summary(m)$r.squared, digits = 3))
  as.character(as.expression(eq))
  
}
lm_eqn(df.t)
r2 <- lm_eqn(df.t)
p1 <-
  p + annotate(
    "text",
    x = 3,
    y = 5,
    label = paste("~italic(r)^2==", r2),
    parse = TRUE,
    color = "gray29",
    size = 7
  )
p1 + labs(
  title = "Gong vs. Gonzalez",
  x = expression(log[2](FC) ~ "Gonzalez"),
  y = expression(log[2](FC) ~ "Gong")
) +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))

dev.off()
dev.off()

# If the gene in the Gong et al paper is also present in the first or full
# term placentas, then keep that gene
Sup4_GongGonzalez <- intersect(Sup4_Gong$Geneid, Gonzalez$Geneid)
Sup4_GongWilson <- intersect(Sup4_Gong$Geneid, Wilson$Geneid)

# Wilson and Gong
Wilson_Sup4 <- subset(Wilson, Wilson$Geneid %in% Sup4_GongWilson)
Gong_Sup4 <-
  subset(Sup4_Gong, Sup4_Gong$Geneid %in% Sup4_GongWilson)

df <- merge(Wilson_Sup4, Gong_Sup4, by = "Geneid")
df.Y <- subset(df, chr == "chrY") #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(1, nrow(df.Y)))
colnames(df.Y)[11] <- "color"
df.X <- subset(df, chr == "chrX") #define Y-chromosomal, blue
df.X <- cbind(df.X, rep(2, nrow(df.X)))
colnames(df.X)[11] <- "color"
df.A <- subset(df, chr != "chrX" & chr != "chrY" ) #define Y-chromosomal, blue
df.A <- cbind(df.A, rep(3, nrow(df.A)))
colnames(df.A)[11] <- "color"
df.t <- rbind(df.Y, df.X, df.A)
df.t$color <- as.factor(df.t$color)


png(
  filename = "FIGURES/Venn_Correlation_DEGs_datasets/log2FC_correlation_Gong_vs_Wilson_Sup4_SDE.png",
  width = 4,
  height = 4,
  units = "in",
  res = 1200
)

p <- ggplot(data = df.t, aes(x = logFC.x, y = logFC.y, color=color)) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "gray29",
    formula = y ~ x,
    size = .5
  ) +
  geom_point(size = 4) + xlim(-2, 8.5) + ylim(-2, 8.5) +
  scale_color_manual(values = c("#66A61E","#7570B3"))#"#D95F02",
p
lm_eqn <- function(df.t) {
  m <- lm(logFC.y ~ logFC.x, df.t)
  
  eq <- list(r2 = format(summary(m)$r.squared, digits = 3))
  as.character(as.expression(eq))
  
}
lm_eqn(df.t)
r2 <- lm_eqn(df.t)
p1 <-
  p + annotate(
    "text",
    x = 3,
    y = 5,
    label = paste("~italic(r)^2==", r2),
    parse = TRUE,
    color = "gray29",
    size = 7
  )
p1 + labs(
  title = "Gong vs. Wilson",
  x = expression(log[2](FC) ~ "Wilson"),
  y = expression(log[2](FC) ~ "Gong")
) +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))

dev.off()
dev.off()

#--- Gonzalez and Gong
Gonzalez_Sup4 <- subset(Gonzalez, Gonzalez$Geneid %in% Sup4_GongGonzalez)
Gong_Sup4 <-
  subset(Sup4_Gong, Sup4_Gong$Geneid %in% Sup4_GongGonzalez)

df <- merge(Gonzalez_Sup4, Gong_Sup4, by = "Geneid")
df.Y <- subset(df, chr == "chrY") #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(1, nrow(df.Y)))
colnames(df.Y)[11] <- "color"
df.X <- subset(df, chr == "chrX") #define Y-chromosomal, blue
df.X <- cbind(df.X, rep(2, nrow(df.X)))
colnames(df.X)[11] <- "color"
df.A <- subset(df, chr != "chrX" & chr != "chrY" ) #define Y-chromosomal, blue
df.A <- cbind(df.A, rep(3, nrow(df.A)))
colnames(df.A)[11] <- "color"
df.t <- rbind(df.Y, df.X, df.A)
df.t$color <- as.factor(df.t$color)


png(
  filename = "FIGURES/Venn_Correlation_DEGs_datasets/log2FC_correlation_Gong_vs_Gonzalez_Sup4_SDE.png",
  width = 4,
  height = 4,
  units = "in",
  res = 1200
)

p <- ggplot(data = df.t, aes(x = logFC.x, y = logFC.y, color=color)) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "gray29",
    formula = y ~ x,
    size = .5
  ) +
  geom_point(size = 4) + xlim(-2, 8.5) + ylim(-2, 8.5) +
  scale_color_manual(values = c("#66A61E","#7570B3"))#"#D95F02",
p
lm_eqn <- function(df.t) {
  m <- lm(logFC.y ~ logFC.x, df.t)
  
  eq <- list(r2 = format(summary(m)$r.squared, digits = 3))
  as.character(as.expression(eq))
  
}
lm_eqn(df.t)
r2 <- lm_eqn(df.t)
p1 <-
  p + annotate(
    "text",
    x = 3,
    y = 7,
    label = paste("~italic(r)^2==", r2),
    parse = TRUE,
    color = "gray29",
    size = 7
  )
p1 + labs(
  title = "Gong vs. Gonzalez",
  x = expression(log[2](FC) ~ "Gonzalez"),
  y = expression(log[2](FC) ~ "Gong")
) +
  theme(plot.title = element_text(size = 12)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))

dev.off()
dev.off()


