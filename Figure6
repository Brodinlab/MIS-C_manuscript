# AutoAb analyses (Figure 6)

# Load IDs
IDs <- read.csv("https://ki.box.com/shared/static/klnb3bnxlbb1uwq94p8i1wvjszgjcqfd.csv", header = T, row.names = 1)
# Load autoAb data (unnormalized signal intensity)
autoAbs <- read.csv("https://ki.box.com/shared/static/m1sudpywenucku4orsfudg3zuc05i7a3.csv", header = T, stringsAsFactors = FALSE)
library(stringr)
autoAbs$uORF <- sapply(str_split(autoAbs$Name, "~"), "[", 3)
autoAbs$uORF <- str_remove(autoAbs$uORF, "uORF:")
autoAbs$duplicates <- paste(autoAbs$accNo, autoAbs$uORF, autoAbs$Block, autoAbs$Row, sep = "_")
autoAbs <- autoAbs[,-c(1:7,9,11:13)] 
# Collect autoAb target info
targetinfo <- autoAbs[,c(1:2,51:52)]
targetinfo <- unique(targetinfo)
# Average duplicated values, final df with 9341 unique autoAbs
autoAbs_avr <- sapply(autoAbs[,-c(1:2,51:52)], tapply, INDEX=autoAbs$duplicates, mean, na.rm=T) 
autoAbs_avr <- as.data.frame(t(autoAbs_avr))
rownames(autoAbs_avr) <- str_remove(rownames(autoAbs_avr), "X")
autoAbs_avr$Row.names <- rownames(autoAbs_avr) 
# Add in sample info to df
autoAbs_sampleinfo <- merge(IDs, autoAbs_avr, by.x="Identification", by.y="Row.names", sort= FALSE)


# Figure 6A
# Df without APS1
autoAbs_noAPS <- autoAbs_sampleinfo[!(autoAbs_sampleinfo$Group %in% "APS1"),]
library(reshape2)
m_autoAbs <- melt(autoAbs_noAPS, id.vars = c("Identification", "Group", "Subjects"))
m_autoAbs$Group <- factor(m_autoAbs$Group, levels = c("Healthy", "CoV2+", "MIS-C", "Kawasaki"))
m_autoAbs <- merge(m_autoAbs, targetinfo, by.x="variable", by.y="duplicates", sort=FALSE)
m_autoAbs$Subjects <- factor(m_autoAbs$Subjects, 
                             levels = c("002917","002865","003073","003080","003079","003236","002984","002963","002973","002911","003222",
                                        "CACTUS 018", "CACTUS 031", "CACTUS 017",
                                        "CACTUS 023", "CACTUS 004", "CACTUS 032",
                                        "Kawasaki_15", "Kawasaki_20", "Kawasaki_16", "Kawasaki_21", "Kawasaki_14", "Kawasaki_19",
                                        "Kawasaki_13", "Kawasaki_17", "Kawasaki_23", "Kawasaki_27", "Kawasaki_31","Kawasaki_24", "Kawasaki_28", "Kawasaki_32",
                                        "Kawasaki_26", "Kawasaki_30", "Kawasaki_25", "Kawasaki_29", "Kawasaki_33", "Kawasaki_3",  "Kawasaki_7",
                                        "Kawasaki_4",  "Kawasaki_8", "Kawasaki_1",  "Kawasaki_6",  "Kawasaki_11", "Kawasaki_5",  "Kawasaki_9" ))

# Violin plots of log value
ggplot(m_autoAbs, aes(y=log(value), x=Subjects)) + 
  geom_violin(aes(fill=Group, colour=Group)) +
  theme_minimal()  + 
  theme(legend.title = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() ) +
  scale_fill_manual(values=c("#2F8AC4", "#E48725", "#A5AA99", "#CD3A8E")) +
  scale_color_manual(values=c("#2F8AC4", "#E48725", "#A5AA99", "#CD3A8E")) +
  labs(x=NULL, y="Log (unnormalized signal)")



# Figure 6B
# CovKD and healthy for GSEA
autoAbs_CovKD <- autoAbs_sampleinfo[autoAbs_sampleinfo$Group %in% c("MIS-C", "Healthy"),]
autoAbs_CovKD$Group <- droplevels(autoAbs_CovKD$Group)
# Calculate min 
df_min <- sapply(autoAbs_CovKD[,-c(1:3)], tapply, INDEX=autoAbs_CovKD$Group, min, na.rm=T)
df_min <- as.data.frame(t(df_min))
colnames(df_min)[1] <- "Healthy_min"
colnames(df_min)[2] <- "MISC_min"
# Calculate max
df_max <- sapply(autoAbs_CovKD[,-c(1:3)], tapply, INDEX=autoAbs_CovKD$Group, max, na.rm=T)
df_max <- as.data.frame(t(df_max))
colnames(df_max)[1] <- "Healthy_max"
colnames(df_max)[2] <- "MISC_max"
# Merge, calculate fold change MIS-C vs Healthy calculate difference MIS-C vs Healthy
df <- merge(df_min, df_max, by="row.names", sort = FALSE)
df$FCminCovKDtomaxH <- df$MISC_min/df$Healthy_max
df <- df[order(df$FCminCovKDtomaxH, decreasing = T),]
df <- merge(df, targetinfo, by.x="Row.names", by.y="duplicates", sort=FALSE)
df <- df[!(duplicated(df$ensembl_ID)),] #remove duplicated ensembl as gse GO uses ensembl identifiers
df$diff <- df$MISC_min-df$Healthy_max #calculate diff between cov kd and healthy
# Ranking based on fold change and difference between MISC and healthy
set.seed(1221)
df$rank2 <- order(order(df$FCminCovKDtomaxH, df$diff, runif(length(df$FCminCovKDtomaxH))))
# GSEA with ranking from fold change and diff
df <- df[order(df$rank2, decreasing = T),]
V2 <- df[,c("ensembl_ID", "rank2")]
library(tibble)
V2 <- deframe(V2)
GGO_V2 <- gseGO(V2, 
                ont="BP",
                keyType = 'ENSEMBL',
                nPerm=1000,
                minGSSize = 3,
                maxGSSize = 500,
                pvalueCutoff = 0.3,
                verbose = T,
                OrgDb='org.Hs.eg.db') # pAdjustMethod =  'BH'
GGO_V2 <- clusterProfiler::simplify(GGO_V2) #removes redundant terms
# GSEA dotplot
dotplot(GGO_V2, orderBy = "GeneRatio") 

# Figure 6C GSEA plot for cardiac ventricle development
gseaplot(GGO_V2, geneSetID = "GO:0003231") 

# Figure 6D 
# Df with GSEA results
gsecc <- GGO_V2[,1:11]
# Genes under cardiac ventricle development GO:0003231
cardiac_ventricle <- gsecc$core_enrichment[gsecc$ID == "GO:0003231"]
cardiac_ventricle <- unlist(str_split(cardiac_ventricle, "/"))
cardiac_ventricle_dupl <- targetinfo$duplicates[targetinfo$ensembl_ID %in% cardiac_ventricle] 
# Violin plots of autoAb targets in cardiac ventricle development GO:0003231
autoAbs_top <- autoAbs_sampleinfo[,c("Identification", "Subjects", "Group", cardiac_ventricle_dupl)]
autoAbs_top <- autoAbs_top[autoAbs_top$Group %in% c("MIS-C", "Healthy", "CoV2+", "Kawasaki"),] #remove APS1 for plotting
library(reshape2)
melt_top <- melt(autoAbs_top, id.vars = c("Identification", "Subjects", "Group"))
melt_top$Group <- factor(melt_top$Group, levels = c("Healthy", "CoV2+", "MIS-C", "Kawasaki"))
melt_top <- merge(melt_top, targetinfo, by.x="variable", by.y="duplicates", sort=FALSE)
melt_top$symbols2 <- paste(melt_top$symbols, melt_top$uORF, sep = "_")
melt_top$symbols2 <- factor(melt_top$symbols2, levels = unique(melt_top$symbols2))
library(ggplot2)
ggplot(melt_top, aes(y=value, x=Group)) + 
  geom_violin(aes(fill=Group, colour=Group)) +
  geom_jitter(shape=16, position=position_jitter(0.2), size=0.8, colour="black") +
  theme_minimal()  + theme(legend.title = element_blank()) +
  scale_fill_manual(values=c("#2F8AC4", "#E48725", "#A5AA99", "#CD3A8E", "white")) +
  scale_color_manual(values=c("#2F8AC4", "#E48725", "#A5AA99", "#CD3A8E", "black")) +
  facet_wrap(~ symbols2, ncol = 7, scales = "free") + labs(x=NULL, y="unnormalized signal") +
  ggtitle("Cardiac ventricle development (GO:0003231)")



