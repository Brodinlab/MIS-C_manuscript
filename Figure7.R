# The Immunology of Multisystem Inflammatory Syndrome in Children with COVID-19
# Figure 7

# Download data from Mendeley Data in Figure_7 folder.
# http://dx.doi.org/10.17632/ds6g796xyg.1


# Preprocessing
# Sample IDs
IDs <- read.csv("Figure7-ID.csv", header = T)

# Quantile-normalized autoAb target dataset
autoAbs <- read.csv("MIS_C.Experiment.QuantileNormailized.200712.csv", header = T, stringsAsFactors = FALSE)
library(stringr)
autoAbs$uORF <- sapply(str_split(autoAbs$Name, "~"), "[", 3)
autoAbs$uORF <- str_remove(autoAbs$uORF, "uORF:")
autoAbs$duplicates <- paste(autoAbs$accNo, autoAbs$uORF, autoAbs$Block, autoAbs$Row, sep = "_")
# AutoAb target info
targetinfo <- autoAbs[,c(7,8,10,71,72)]
targetinfo <- unique(targetinfo)
# Average duplicated values, final df with 9341 unique autoAbs
autoAbs_avr <- as.data.frame(sapply(autoAbs[,14:70], tapply, INDEX=autoAbs$duplicates, mean, na.rm=T))
colnames(autoAbs_avr) <- str_remove(colnames(autoAbs_avr), "X")
autoAbs_avr <- as.data.frame(t(autoAbs_avr))
rm(autoAbs)
autoAbs_avr <- autoAbs_avr[!(rownames(autoAbs_avr) == "90021"),] #Remove duplicated sample MISC 16 (90021), which is a duplicate of MISC 14
# Add in sample info to df
autoAbs_sampleinfo <- merge(IDs, autoAbs_avr, by.x="Identification", by.y="row.names", sort= FALSE)
rownames(autoAbs_sampleinfo) <- autoAbs_sampleinfo$Identification
# Melt
library(reshape2)
m_autoAbs <- reshape2::melt(autoAbs_sampleinfo, id.vars = c("Identification", "Subjects", "Subjects2", "Group",  "Group2", "Batch"))
m_autoAbs$Group <- factor(m_autoAbs$Group, levels = c("Healthy", "CoV2+", "MIS-C", "Kawasaki"))
m_autoAbs <- merge(m_autoAbs, targetinfo, by.x="variable", by.y="duplicates", sort=FALSE)
m_autoAbs$symbols2 <- paste(m_autoAbs$symbols, m_autoAbs$uORF, sep = "_")
mylevelsS1 <- c("2917","2865","3073","3080","3079","3236","2984","2963","2973","2911","3222",
                "K1", "K2",
                "CACTUS 018", "CACTUS 031", "CACTUS 017",
                "CACTUS 023", "CACTUS 004", "CACTUS 032",
                "MIS-C 2", "MIS-C 3", "MIS-C 6", "MIS-C 8", "MIS-C 10", "MIS-C 12", "MIS-C 14",
                "CV5", "CV 9",
                "Kawasaki_15", "Kawasaki_20", "Kawasaki_16", "Kawasaki_21", "Kawasaki_14", "Kawasaki_19",
                "Kawasaki_13", "Kawasaki_17", "Kawasaki_23", "Kawasaki_27", "Kawasaki_31","Kawasaki_24", "Kawasaki_28", "Kawasaki_32",
                "Kawasaki_26", "Kawasaki_30", "Kawasaki_25", "Kawasaki_29", "Kawasaki_33", "Kawasaki_3", "Kawasaki_7",
                "Kawasaki_4",  "Kawasaki_8", "Kawasaki_1", "Kawasaki_6", "Kawasaki_11", "Kawasaki_5", "Kawasaki_9")
m_autoAbs$Subjects <- factor(m_autoAbs$Subjects, levels = mylevelsS1)

# Figure 7A
library(ggplot2)
ggplot(m_autoAbs, aes(y=log1p(value), x=Subjects)) + 
  geom_violin(aes(fill=Group, colour=Group)) +
  theme_minimal()  + 
  theme(legend.title = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() ) +
  scale_fill_manual(values=c("#2F8AC4", "#E48725", "#A5AA99", "#CD3A8E")) +
  scale_color_manual(values=c("#2F8AC4", "#E48725", "#A5AA99", "#CD3A8E")) +
  labs(x=NULL, y="Log (1 + Normalized signal)")


# Data analysis
# Deseq: MISC vs groups
txi <- merge(IDs, autoAbs_avr, by.x="Identification", by.y="row.names")
rownames(txi) <- txi$Identification
meta <- txi[,1:6]
meta[,c(2,4,5,6)] <- sapply(meta[,c(2,4,5,6)], function(i) str_replace(i, "\\-", "_"))
meta[,c(2,4,5,6)] <- sapply(meta[,c(2,4,5,6)], function(i) str_replace(i, " ", "_"))
meta[,c(2,4,5,6)] <- sapply(meta[,c(2,4,5,6)], function(i) str_replace(i, "\\+", ""))
txi <- txi[,-(1:6)]
txi <- as.data.frame(t(txi))
#identical(rownames(meta), colnames(txi)) #TRUE
txi[sapply(txi, function(i) i < 0)] <- 0
txi[,1:56] <- sapply(txi[,1:56], as.integer)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = txi, colData = meta, design = ~ Batch + Group)
dds$Batch <- factor(dds$Batch, levels = c("A","B"))
dds$Group <- factor(dds$Group, levels = c("MIS_C","Healthy","CoV2", "Kawasaki")) #we want to compare MISC with others
dds <- estimateSizeFactors(dds)
dds1 <- DESeq(dds, fitType='local')
# Results
# MISC vs Healthy
Healthy_MISC <- results(dds1, name="Group_Healthy_vs_MIS_C") 
Healthy_MISC_1 <- as.data.frame(Healthy_MISC)
Healthy_MISC_1 <- merge(Healthy_MISC_1, targetinfo, by.x="row.names", by.y="duplicates", sort=FALSE)
# MISC vs CoV2
CoV_MISC <- results(dds1, name="Group_CoV2_vs_MIS_C") 
CoV_MISC_1 <- as.data.frame(CoV_MISC)
CoV_MISC_1 <- merge(CoV_MISC_1, targetinfo, by.x="row.names", by.y="duplicates", sort=FALSE)
# Kawasaki vs CoV2
Kawasaki_MISC <- results(dds1, name="Group_Kawasaki_vs_MIS_C") 
Kawasaki_MISC_1 <- as.data.frame(Kawasaki_MISC)
Kawasaki_MISC_1 <- merge(Kawasaki_MISC_1, targetinfo, by.x="row.names", by.y="duplicates", sort=FALSE)

# GSEA 
library(tibble)
library(clusterProfiler)
library(org.Hs.eg.db)
# GSEA MISC vs Healthy
Healthy_MISC_1$rank <- -Healthy_MISC_1$log2FoldChange
Healthy_MISC_1 <- Healthy_MISC_1[order(Healthy_MISC_1$rank, decreasing = TRUE),] 
Healthy_MISC_1 <- Healthy_MISC_1[!duplicated(Healthy_MISC_1$ensembl_ID),]
rownames(Healthy_MISC_1) <- NULL
Healthy_MISC_2 <- Healthy_MISC_1[,c("ensembl_ID", "rank")]
Healthy_MISC_2 <- deframe(Healthy_MISC_2)
set.seed(1234)
GGO_HvsM <- gseGO(Healthy_MISC_2, 
                  ont="BP",
                  keyType = 'ENSEMBL',
                  nPerm=1000,
                  minGSSize = 60, 
                  maxGSSize = 500, 
                  pvalueCutoff = 0.05,
                  verbose = T,
                  OrgDb='org.Hs.eg.db',
                  seed = TRUE) 
GGO_HvsM <- clusterProfiler::simplify(GGO_HvsM)
GGO_HvsM_df <- GGO_HvsM[,1:11]
GGO_HvsM_df <- GGO_HvsM_df[order(GGO_HvsM_df$NES, decreasing = TRUE),]

# GSEA MISC vs CoV2+
CoV_MISC_1$rank <- -CoV_MISC_1$log2FoldChange
CoV_MISC_1 <- CoV_MISC_1[order(CoV_MISC_1$rank, decreasing = TRUE),] 
CoV_MISC_1 <- CoV_MISC_1[!duplicated(CoV_MISC_1$ensembl_ID),]
rownames(CoV_MISC_1) <- NULL
CoV_MISC_2 <- CoV_MISC_1[,c("ensembl_ID", "rank")]
CoV_MISC_2 <- deframe(CoV_MISC_2)
set.seed(1234)
GGO_CvsM <- gseGO(CoV_MISC_2, 
                  ont="BP",
                  keyType = 'ENSEMBL',
                  nPerm=1000,
                  minGSSize = 60, 
                  maxGSSize = 500, 
                  pvalueCutoff = 0.05,
                  verbose = T,
                  OrgDb='org.Hs.eg.db',
                  seed = TRUE) 
GGO_CvsM <- clusterProfiler::simplify(GGO_CvsM)
GGO_CvsM_df <- GGO_CvsM[,1:11]
GGO_CvsM_df <- GGO_CvsM_df[order(GGO_CvsM_df$NES, decreasing = TRUE),]

# GSEA MISC vs Kawasaki
Kawasaki_MISC_1$rank <- -Kawasaki_MISC_1$log2FoldChange
Kawasaki_MISC_1 <- Kawasaki_MISC_1[order(Kawasaki_MISC_1$rank, decreasing = TRUE),] 
Kawasaki_MISC_1 <- Kawasaki_MISC_1[!duplicated(Kawasaki_MISC_1$ensembl_ID),]
rownames(Kawasaki_MISC_1) <- NULL
Kawasaki_MISC_2 <- Kawasaki_MISC_1[,c("ensembl_ID", "rank")]
Kawasaki_MISC_2 <- deframe(Kawasaki_MISC_2)
set.seed(1234)
GGO_KvsM <- gseGO(Kawasaki_MISC_2, 
                  ont="BP",
                  keyType = 'ENSEMBL',
                  nPerm=1000,
                  minGSSize = 60, 
                  maxGSSize = 500, 
                  pvalueCutoff = 0.05,
                  verbose = T,
                  OrgDb='org.Hs.eg.db',
                  seed = TRUE) 
GGO_KvsM <- clusterProfiler::simplify(GGO_KvsM)
GGO_KvsM_df <- GGO_KvsM[,1:11]
GGO_KvsM_df <- GGO_KvsM_df[order(GGO_KvsM_df$NES, decreasing = TRUE),]

# VennDiagram
library(VennDiagram)
y <- list()
y$MISCvsHealthy <- paste(GGO_HvsM_df$ID, GGO_HvsM_df$Description, sep = " ")
y$MISCvsCoV2 <- paste(GGO_CvsM_df$ID, GGO_CvsM_df$Description, sep = " ")
y$MISCvsKawasaki <- paste(GGO_KvsM_df$ID, GGO_KvsM_df$Description, sep = " ")
v0 <-venn.diagram(y, lwd = 3, col = c("black", "black", "black"), fill = c("#2F8AC4", "#E48725", "#CD3A8E"), alpha = 0.5, filename = NULL, cex=1.5)
# Figure 7B
grid.newpage()
grid.draw(v0)

# Intersection of GO terms between analyses
int1 <- intersect(y$MISCvsHealthy, y$MISCvsCoV2) 
int2 <- intersect(y$MISCvsHealthy, y$MISCvsKawasaki) 
int3 <- intersect(y$MISCvsCoV2, y$MISCvsKawasaki) 
int4 <- intersect(int1, int2) # 26 common to all comparisons
# Plot enriched go terms common to all MISC comparisons 
#library(UpSetR)
library(ComplexHeatmap)
z <- list_to_matrix(y)
z <- as.data.frame(z)
zz <- z[rownames(z) %in% int4,]
zz$rownames <- rownames(zz)
zz <- reshape2::melt(zz, id.vars="rownames")

# Figure 7C
ggplot(na.omit(zz), aes(x=variable, y=rownames)) + 
  geom_point() + ylab(NULL) + xlab(NULL) + 
  theme_minimal() + theme(legend.position = "none") + theme(axis.text.y = element_text(size = 5))

# Figure 7D
# GO:0007507 heart develpment MIS-C vs CoV2+
gseaplot(GGO_CvsM, "GO:0007507")


# Figure 7E
ENG <-Healthy_MISC_1$Row.names[Healthy_MISC_1$symbols %in% "ENG"]
autoAbs_top <- autoAbs_sampleinfo[,c("Identification", "Subjects", "Group", ENG)]
autoAbs_top <- autoAbs_top[autoAbs_top$Group %in% c("MIS-C", "Healthy", "CoV2+", "Kawasaki"),] 
melt_top <- reshape2::melt(autoAbs_top, id.vars = c("Identification", "Subjects", "Group"))
melt_top$Group <- factor(melt_top$Group, levels = c("Healthy", "CoV2+", "MIS-C", "Kawasaki"))
melt_top <- merge(melt_top, targetinfo, by.x="variable", by.y="duplicates", sort=FALSE)
melt_top$symbols <- factor(melt_top$symbols, levels = unique(melt_top$symbols))
ggplot(melt_top, aes(y=value, x=Group)) + 
  geom_violin(aes(fill=Group, colour=Group)) +
  geom_jitter(shape=16, position=position_jitter(0.2), size=0.8, colour="black") +
  theme_minimal()  + theme(legend.title = element_blank()) +
  scale_fill_manual(values=c("#2F8AC4", "#E48725", "#A5AA99", "#CD3A8E")) +
  scale_color_manual(values=c("#2F8AC4", "#E48725", "#A5AA99", "#CD3A8E")) +
  facet_wrap(~ symbols, ncol = 7, scales = "free") + labs(x=NULL, y="Normalized signal") 
pairwise.wilcox.test(x = melt_top$value, g = melt_top$Group, p.adjust.method = "fdr")


# Figure 7F
ELISA <- read.csv("ELISA.csv", stringsAsFactors = FALSE, row.names = 1)
ELISA$Group <- factor(ELISA$Group, levels = c("Healthy", "CoV2+", "MIS-C", "Kawasaki"))
ggplot(ELISA, aes(x=Group, y=log(pg.mL))) + theme_minimal() +
  geom_boxplot(aes(fill=Group)) +
  geom_jitter(shape=16, position=position_jitter(0.2), size=0.8, colour="black") +
  scale_fill_manual(values=c("#2F8AC4", "#E48725", "#A5AA99", "#CD3A8E")) +
  theme(axis.text.x = element_text(size = 5), legend.title = element_blank()) + labs(x=NULL, y="Log(Endoglin pg/mL)") 
pairwise.wilcox.test(ELISA$pg.mL, ELISA$Group, p.adjust.method = "fdr")


# Figure 7G
serine <- GGO_CvsM_df$core_enrichment[GGO_CvsM_df$ID == "GO:0018105"]
serine <- unlist(str_split(serine, "/"))
serine_dupl <- CoV_MISC_1$Row.names[CoV_MISC_1$ensembl_ID %in% serine]
autoAbs_top <- autoAbs_sampleinfo[,c("Identification", "Subjects", "Group", serine_dupl)]
autoAbs_top <- autoAbs_top[autoAbs_top$Group %in% c("MIS-C", "Healthy", "CoV2+", "Kawasaki"),] 
melt_top <- reshape2::melt(autoAbs_top, id.vars = c("Identification", "Subjects", "Group"))
melt_top$Group <- factor(melt_top$Group, levels = c("Healthy", "CoV2+", "MIS-C", "Kawasaki"))
melt_top <- merge(melt_top, targetinfo, by.x="variable", by.y="duplicates", sort=FALSE)
melt_top <- melt_top[melt_top$symbols %in% c("MAP2K2", "CSNK1A1", "CSNK2A1", "CSNK1E"),]
ggplot(melt_top, aes(y=value, x=Group)) + 
  geom_violin(aes(fill=Group, colour=Group)) +
  geom_jitter(shape=16, position=position_jitter(0.2), size=0.8, colour="black") +
  theme_minimal()  + theme(legend.title = element_blank()) +
  scale_fill_manual(values=c("#2F8AC4", "#E48725", "#A5AA99", "#CD3A8E")) +
  scale_color_manual(values=c("#2F8AC4", "#E48725", "#A5AA99", "#CD3A8E")) +
  facet_wrap(~ symbols, ncol = 1, scales = "free") + labs(x=NULL, y="Normalized signal") +
  ggtitle("GO:0018105 Peptidyl-serine phosphorylation")

# Figure 7H
# Volcano plot Kawasaki vs MISC
# Kawasaki vs CoV2
Kawasaki_MISC <- results(dds1, name="Group_Kawasaki_vs_MIS_C") 
Kawasaki_MISC_1 <- as.data.frame(Kawasaki_MISC)
Kawasaki_MISC_1 <- merge(Kawasaki_MISC_1, targetinfo, by.x="row.names", by.y="duplicates", sort=FALSE)
Kawasaki_MISC_1$Sig <- ifelse(Kawasaki_MISC_1$padj < 0.05, Kawasaki_MISC_1$Sig <- "pAdj < 0.05", Kawasaki_MISC_1$Sig <- "ns")
library(ggrepel)
ggplot(subset(Kawasaki_MISC_1, !is.na(Kawasaki_MISC_1$padj)), aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = Sig)) +
  scale_color_manual(values = c("grey", "#CD3A8E")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  labs(y="-Log10(pValue)", x="Log2 (Kawasaki/MIS-C)") +
  xlim(c(-10,15)) + ylim(c(0,20)) +
  geom_text_repel(
    data = subset(Kawasaki_MISC_1, padj < 0.001), 
    aes(label = symbols),
    size = 2,
    segment.size = 0.05,
    segment.colour = "transparent") 


sessionInfo()
#R version 3.6.2 (2019-12-12)
#Platform: x86_64-apple-darwin15.6.0 (64-bit)
#Running under: macOS Catalina 10.15.6

#attached base packages:
#[1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base   

#other attached packages:
#[1] ggrepel_0.8.2               ComplexHeatmap_2.0.0        VennDiagram_1.6.20          futile.logger_1.4.3        
#[5] org.Hs.eg.db_3.8.2          AnnotationDbi_1.46.1        clusterProfiler_3.12.0      tibble_3.0.3               
#[9] DESeq2_1.24.0               SummarizedExperiment_1.14.1 DelayedArray_0.10.0         BiocParallel_1.18.1        
#[13] matrixStats_0.56.0          Biobase_2.44.0              GenomicRanges_1.36.1        GenomeInfoDb_1.20.0        
#[17] IRanges_2.18.3              S4Vectors_0.22.1            BiocGenerics_0.30.0         ggplot2_3.3.2              
#[21] reshape2_1.4.4              stringr_1.4.0





