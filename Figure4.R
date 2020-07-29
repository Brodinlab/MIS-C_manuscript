# The Immunology of Multisystem Inflammatory Syndrome in Children with COVID-19
# Figure 4

# Download olink data from:
# https://ki.box.com/s/7qwdagyeqo7yhs6502aaxxb06oc5otlp

# Plasma protein data preprocessing
# Read Olink NPX data
library(OlinkAnalyze)
May <- read_NPX(filename = "OlinkMay1.xlsx") #  NPX values under LOD already substituted by LOD
July <- read_NPX(filename = "OlinkJuly1.xlsx")
July$NPX[July$NPX <= July$LOD] <- July$LOD[July$NPX <= July$LOD] # substitute NPX values under LOD for LOD
# Correct NPX values for samples that were diluted
July$unnormalized <- 2^(July$NPX)
dilutedby2 <- c("CACTUS 004", "CACTUS 031", "CACTUS 032", "CACTUS 018", "CACTUS 023", "CACTUS 024", "CACTUS 017")
July$dilutionfactor <- ifelse(July$SampleID %in% dilutedby2, July$dilutionfactor <- 2, July$dilutionfactor <- 1 )
July$unnormalized2 <- (July$unnormalized)*(July$dilutionfactor)
July$NPX2 <- log2(July$unnormalized2)
newJuly <- July
newJuly$NPX <- newJuly$NPX2  
newJuly <- newJuly[,-c(14:17)]
rm(July, dilutedby2)
# bridge normalization (using bridge samples common in both runs)
bridge_samples <- intersect(x = newJuly$SampleID, y = May$SampleID)
bridge_normalized_data <- olink_normalization(df1 = May, df2 = newJuly, overlapping_samples_df1 = bridge_samples)
# cast df
olink <- bridge_normalized_data[,c(1,5,12)]
library(reshape2)
olink <- dcast(olink, SampleID ~ factor(Assay, levels = unique(olink$Assay)), fun.aggregate=function(i) mean(i, na.rm=TRUE))
rownames(olink) <- olink$SampleID
olink <- olink[,-1]
# get LOD; since LOD varies between panels, collect individual LOD per protein then average them out
LOD <- bridge_normalized_data[,c(5,8,11)]
LOD$dupl <- paste(LOD$Assay, LOD$Panel_Version, sep = "__")
LOD <- as.data.frame(LOD[!duplicated(LOD$dupl),])
LOD <- LOD[,-c(2,4)]
library(dplyr)
LOD$Assay <- factor(LOD$Assay, levels = unique(LOD$Assay))
LOD <- as.data.frame(LOD %>% group_by(Assay) %>% summarise(mean=mean(LOD)))
LOD$Assay <- as.character(LOD$Assay)
rownames(LOD) <- LOD$Assay
identical(LOD$Assay,colnames(olink)) #TRUE
#add LOD to olink df
olink <- as.data.frame(t(olink))
olink <- merge(olink, LOD[,2, drop=FALSE], by="row.names", all = TRUE)
rownames(olink) <- olink$Row.names
olink$Row.names <- NULL
colnames(olink)[106] <- "LOD"
# Detect number of samples with limit of detection (LOD)
num_LOD <- as.data.frame(apply(olink, 1, function(x) length(which(x[1:105] <= x[106])) ))
olink <- merge(olink, num_LOD, by="row.names")
rownames(olink) <- olink$Row.names
olink <- olink[,-1]
colnames(olink)[107] <- "n_LOD"
olink <- olink[,-106] #remove LOD
olink <- as.data.frame(t(olink))
# Remove proteins from analysis where more than 30% of samples have LOD 
# 30% of 105 is 31.5
olink <- olink[,-which(olink[106,] > 31.5)] # 133 proteins
olink <- olink[-106,] #remove lod_n info
# Read ID info
ID <- read.csv("https://ki.box.com/shared/static/vcrn9uaud6b5l7uk5rvgjgkf49mvehcv.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
# Merge olink with ID data
olink <- merge(ID[,c(1,2)], olink, by.x="Sample", by.y="row.names", sort=FALSE)
MyBatch <- olink[,1:2]
rownames(olink) <- olink$Sample
olink <- olink[-c(1:2)]
olink <- t(olink)
# Batchcorrect
library(limma)
olinkX <- as.data.frame(t(removeBatchEffect(x = olink, batch = MyBatch$Batch2)))
# Remove healthy adults from olink df
olinkX <- subset(olinkX, !(grepl("HC", rownames(olinkX) )) ) #102 subjects total
# Remove multiple samples from same subject
# MISC6 and MISC5 are the same subject/ same day / pretreatment
# MISC11 and MISC12 are the same subject/ same day / pretreatment
# MISC15 and MISC1 are the same subject/ same day / posttreatment
olinkX <- olinkX[!(rownames(olinkX) %in% c("MIS-C 6", "MIS-C 11", "MIS-C 15")),] #99 subjects total
# Remove samples from treated MIS-C subjects
post <- ID$Sample[ID$Treatment == "Post"]
rownames(olinkX)[(rownames(olinkX) %in% post)]
olinkX <- olinkX[!(rownames(olinkX) %in% post),] #92 total
rm(bridge_samples, bridge_normalized_data, LOD, May, MyBatch, newJuly, num_LOD, olink, post)

# Data analysis
# PCA
pca <- prcomp(na.omit(olinkX), scale. = T)
pca_x <- as.data.frame(pca$x)
pca_x <- merge(ID, pca_x, by.x="Sample", by.y="row.names", sort=FALSE)
pca_x$Group2 <- factor(pca_x$Group2, levels = c("Healthy kids", "CoV2+", "MIS-C", "Kawasaki"))
# Scree plot
library(factoextra)
scree <- fviz_eig(pca)
screedata <- scree$data
sum(screedata$eig[1:5]) #Principal comp 1:5 explain 58.3153 % of variation in data

# Figure 4A
# PCA pairs
# Plotting function, color by group
color_by_group <- function(pca_x, xaxis, yaxis) {
  p <- ggplot(pca_x, aes(y = .data[[yaxis]], x = .data[[xaxis]]))
  p + geom_point(aes(color=Group2)) + theme_bw() +
    theme(legend.title = element_blank(), legend.position="none") + 
    scale_color_manual(values=c("#2F8AC4", "#E48725", "#A5AA99", "#CD3A8E")) + 
    labs(x=xaxis, y=yaxis) }

# PCA pairs PC1:PC5
PCA12 <- color_by_group(pca_x, "PC1", "PC2") 
PCA13 <- color_by_group(pca_x, "PC1", "PC3") 
PCA23 <- color_by_group(pca_x, "PC2", "PC3") 
PCA14 <- color_by_group(pca_x, "PC1", "PC4")
PCA24 <- color_by_group(pca_x, "PC2", "PC4")
PCA34 <- color_by_group(pca_x, "PC3", "PC4")
PCA15 <- color_by_group(pca_x, "PC1", "PC5")
PCA25 <- color_by_group(pca_x, "PC2", "PC5")
PCA35 <- color_by_group(pca_x, "PC3", "PC5")
PCA45 <- color_by_group(pca_x, "PC4", "PC5")
library(cowplot)
legend <- get_legend(ggplot(pca_x, aes(x=PC1, y=PC2)) + geom_point(aes(color=Group2)) + theme_bw() +
                       theme(legend.title = element_blank()) + 
                       scale_color_manual(values=c("#2F8AC4", "#E48725", "#A5AA99", "#CD3A8E")) + theme(legend.box.margin = margin(0, 0, 0, 12)) ) # create some space to the left of the legend
plot_grid(PCA12, NULL, NULL, legend, 
          PCA13, PCA23, NULL, NULL, 
          PCA14, PCA24, PCA34, NULL, 
          PCA15, PCA25, PCA35, PCA45,
          ncol = 4)

# Figure 4B
# Top 25 loadings of PC2
library(factoextra)
load_PC2 <- fviz_pca_var(pca,
                         axes = c(2,2),
                         col.var = "contrib", 
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel = TRUE,
                         labelsize=2,
                         select.var = list(contrib = 25) ) 
load_PC2data <- as.data.frame(load_PC2$data)
load_PC2data <- load_PC2data[order(load_PC2data$x, decreasing = T),]
load_PC2data$name <- factor(load_PC2data$name, levels = load_PC2data$name)
ggplot(load_PC2data, aes(x=name, y=x)) +
  geom_segment( aes(x=name, xend=name, y=0, yend=x), color="grey") +
  geom_point( color="#A5AA99", size=4) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank() ) +
  xlab("") +
  ylab("Contribution to PC2")


# Fig 4C, D, E
# Violin plots
olinkX_ID <- merge(ID, olinkX, by.x="Sample", by.y="row.names", sort=FALSE)
melt_olink <- reshape2::melt(olinkX_ID, id.vars = c("Sample", "Group2", "Batch2", "Treatment"))
melt_olink$Group2 <- factor(melt_olink$Group2, levels = c("Healthy kids", "CoV2+","MIS-C", "Kawasaki"))
Fig4C <- c("IL6", "IL-17A", "CXCL10")
Fig4D <- c("SCF", "TWEAK", "ADA")
Fig4E <- c("DCBLD2", "EDAR", "CLEC4C")
Fig4CDE <- c(Fig4C, Fig4D, Fig4E)
melt_olink <- subset(melt_olink, variable %in% Fig4CDE)
melt_olink$variable <- factor(melt_olink$variable, levels = Fig4CDE)
ggplot(melt_olink, aes(x=Group2, y=value)) + theme_minimal() + 
  geom_violin(aes(fill=Group2, colour=Group2)) +
  facet_wrap(~ variable, scales = "free", ncol = 3) + 
  geom_jitter(shape=16, position=position_jitter(0.2), size=0.8, colour="black") +
  scale_fill_manual(values=c("#2F8AC4", "#E48725", "#A5AA99", "#CD3A8E")) +
  scale_color_manual(values=c("#2F8AC4", "#E48725", "#A5AA99", "#CD3A8E")) +
  theme(axis.text.x = element_text(size = 5), legend.title = element_blank()) + labs(x=NULL, y="NPX") 
# pvalues
olink_kd_MISC <- olinkX_ID[olinkX_ID$Group2 %in% c("MIS-C", "Kawasaki"), ]
olink_kd_MISC <- olink_kd_MISC[,c("Sample","Group2",Fig4CDE)]
set.seed(1234)
diffs <- lapply(olink_kd_MISC[,-(1:2)], function(i) wilcox.test(i ~ olink_kd_MISC$Group2)$p.value )
diffs <- unlist(diffs)
diffs


# Supplementary Figure 3
# Supplementary Figure 3A
melt_olink <- reshape2::melt(olinkX_ID, id.vars = c("Sample", "Group2", "Batch2", "Treatment"))
melt_olink$Group2 <- factor(melt_olink$Group2, levels = c("Healthy kids", "CoV2+","MIS-C", "Kawasaki"))
MMPs <- c("MMP-10", "MMP-1")
melt_olink <- subset(melt_olink, variable %in% MMPs)
ggplot(melt_olink, aes(x=Group2, y=value)) + theme_minimal() + 
  geom_violin(aes(fill=Group2, colour=Group2)) +
  facet_wrap(~ variable, scales = "free", ncol = 3) + 
  geom_jitter(shape=16, position=position_jitter(0.2), size=0.8, colour="black") +
  scale_fill_manual(values=c("#2F8AC4", "#E48725", "#A5AA99", "#CD3A8E")) +
  scale_color_manual(values=c("#2F8AC4", "#E48725", "#A5AA99", "#CD3A8E")) +
  theme(axis.text.x = element_text(size = 5), legend.title = element_blank()) + labs(x=NULL, y="NPX") 
# pvalues
olink_kd_MISC <- olinkX_ID[olinkX_ID$Group2 %in% c("MIS-C", "Kawasaki"), ]
olink_kd_MISC <- olink_kd_MISC[,c("Sample","Group2",MMPs)]
set.seed(1234)
diffs <- lapply(olink_kd_MISC[,-(1:2)], function(i) wilcox.test(i ~ olink_kd_MISC$Group2)$p.value )
diffs <- unlist(diffs)
diffs

# Supplementary Figure 3B
olink_kd_CoV <- olinkX_ID[olinkX_ID$Group2 %in% c("MIS-C", "CoV2+"), ]
cdiffs <- lapply(olink_kd_CoV[,-(1:4)], function(i) wilcox.test(i ~ olink_kd_CoV$Group2)$p.value )
cdiffs <- unlist(cdiffs)
cdiffs <- cdiffs[cdiffs < 0.05]
cdiffs
covdiffs <- names(cdiffs)
melt_olink <- reshape2::melt(olinkX_ID, id.vars = c("Sample", "Group2", "Batch2", "Treatment"))
melt_olink$Group2 <- factor(melt_olink$Group2, levels = c("Healthy kids", "CoV2+","MIS-C", "Kawasaki"))
melt_olink <- subset(melt_olink, variable %in% covdiffs)
ggplot(melt_olink, aes(x=Group2, y=value)) + theme_minimal() + 
  geom_violin(aes(fill=Group2, colour=Group2)) +
  facet_wrap(~ variable, scales = "free", ncol = 5) + 
  geom_jitter(shape=16, position=position_jitter(0.2), size=0.8, colour="black") +
  scale_fill_manual(values=c("#2F8AC4", "#E48725", "#A5AA99", "#CD3A8E")) +
  scale_color_manual(values=c("#2F8AC4", "#E48725", "#A5AA99", "#CD3A8E")) +
  theme(axis.text.x = element_text(size = 5), legend.title = element_blank()) + labs(x=NULL, y="NPX") 



sessionInfo()
#R version 3.6.2 (2019-12-12)
#Platform: x86_64-apple-darwin15.6.0 (64-bit)
#Running under: macOS Catalina 10.15.6

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] cowplot_1.0.0      factoextra_1.0.7   limma_3.40.6       reshape2_1.4.4     OlinkAnalyze_1.0.2 readxl_1.3.1       forcats_0.5.0     
#[8] dplyr_1.0.0        purrr_0.3.4        readr_1.3.1        tidyr_1.1.0        tibble_3.0.2       tidyverse_1.3.0    stringr_1.4.0     
#[15] openxlsx_4.1.5     lmerTest_3.1-2     lme4_1.1-23        Matrix_1.2-18      ggrepel_0.8.2      ggfortify_0.4.10   ggplot2_3.3.2     
#[22] emmeans_1.4.8      car_3.0-8          carData_3.0-4      broom_0.7.0       


