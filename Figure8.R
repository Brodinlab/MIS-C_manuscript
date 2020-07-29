# The Immunology of Multisystem Inflammatory Syndrome in Children with COVID-19
# Figure 8

# Preprocessing
# Sample IDs
IDs <- read.csv("https://ki.box.com/shared/static/uq3bns63ebakzqbtzuduqivcyalwiv2k.csv", header = T)
# Quantile-normalized autoAb target dataset
autoAbs <- read.csv("https://ki.box.com/shared/static/z4a4torxyevlx5a00ksp833bi87lofkc.csv", header = T, stringsAsFactors = FALSE)
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
dds$Group <- factor(dds$Group, levels = c("MIS_C","Healthy","CoV2", "Kawasaki")) #MISC vs others
dds <- estimateSizeFactors(dds)
dds1 <- DESeq(dds, fitType='local')

# GSEA 
library(tibble)
library(clusterProfiler)
library(org.Hs.eg.db)
# MISC vs CoV2
CoV_MISC <- results(dds1, name="Group_CoV2_vs_MIS_C") 
CoV_MISC_1 <- as.data.frame(CoV_MISC)
CoV_MISC_1 <- merge(CoV_MISC_1, targetinfo, by.x="row.names", by.y="duplicates", sort=FALSE)
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

# Figure 8A
# GO:0018105	peptidyl-serine phosphorylation
gseaplot(GGO_CvsM, "GO:0018105")

# Figure 8B
# Violin plots heart development
# Add in sample info to df
autoAbs_sampleinfo <- merge(IDs, autoAbs_avr, by.x="Identification", by.y="row.names", sort= FALSE)
rownames(autoAbs_sampleinfo) <- autoAbs_sampleinfo$Identification
serine <- GGO_CvsM_df$core_enrichment[GGO_CvsM_df$ID == "GO:0018105"]
serine <- unlist(str_split(serine, "/"))
serine_dupl <- CoV_MISC_1$Row.names[CoV_MISC_1$ensembl_ID %in% serine]
autoAbs_top <- autoAbs_sampleinfo[,c("Identification", "Subjects", "Group", serine_dupl)]
autoAbs_top <- autoAbs_top[autoAbs_top$Group %in% c("MIS-C", "Healthy", "CoV2+", "Kawasaki"),] 
melt_top <- reshape2::melt(autoAbs_top, id.vars = c("Identification", "Subjects", "Group"))
melt_top$Group <- factor(melt_top$Group, levels = c("Healthy", "CoV2+", "MIS-C", "Kawasaki"))
melt_top <- merge(melt_top, targetinfo, by.x="variable", by.y="duplicates", sort=FALSE)
melt_top <- melt_top[melt_top$symbols %in% c("MAP2K2", "CSNK1A1", "CSNK2A1", "CSNK1E"),]
library(ggplot2)
ggplot(melt_top, aes(y=value, x=Group)) + 
  geom_violin(aes(fill=Group, colour=Group)) +
  geom_jitter(shape=16, position=position_jitter(0.2), size=0.8, colour="black") +
  theme_minimal()  + theme(legend.title = element_blank()) +
  scale_fill_manual(values=c("#2F8AC4", "#E48725", "#A5AA99", "#CD3A8E")) +
  scale_color_manual(values=c("#2F8AC4", "#E48725", "#A5AA99", "#CD3A8E")) +
  facet_wrap(~ symbols, ncol = 1, scales = "free") + labs(x=NULL, y="Normalized signal") +
  ggtitle("GO:0018105 Peptidyl-serine phosphorylation")


# Figure 8 C
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
#[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] ggrepel_0.8.2               ggplot2_3.3.2               org.Hs.eg.db_3.8.2          AnnotationDbi_1.46.1       
#[5] clusterProfiler_3.12.0      tibble_3.0.2                DESeq2_1.24.0               SummarizedExperiment_1.14.1
#[9] DelayedArray_0.10.0         BiocParallel_1.18.1         matrixStats_0.56.0          Biobase_2.44.0             
#[13] GenomicRanges_1.36.1        GenomeInfoDb_1.20.0         IRanges_2.18.3              S4Vectors_0.22.1           
#[17] BiocGenerics_0.30.0         stringr_1.4.0




