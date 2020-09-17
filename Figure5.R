# The Immunology of Multisystem Inflammatory Syndrome in Children with COVID-19
# Figure 5

# Download data from Mendeley Data in Figure_5 folder.
# http://dx.doi.org/10.17632/ds6g796xyg.1

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
ID <- read.csv("Figure5-ID.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
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
olinkX_ID <- merge(ID, olinkX, by.x="Sample", by.y="row.names", sort=FALSE)
rm(bridge_samples, bridge_normalized_data, LOD, May, MyBatch, newJuly, num_LOD, olink)

# Data analysis of paired samples
# MISC untreated vs treated
unt_tr <- olinkX_ID[olinkX_ID$Group2 %in% "MIS-C", ]
MISCpairs <- c("MIS-C 032", "MIS-C 004", "MIS-C 023", "MIS-C 002",  "MIS-C 005",  "MIS-C 008",  "MIS-C 014" )
un_tr_pairs <- unt_tr[unt_tr$Subject %in% MISCpairs,] 
rownames(un_tr_pairs) <- NULL
un_tr_pairs <- un_tr_pairs[c(6,1,3,5,4,2,7,12,8,9,10,13,14,11),] #order samples
diffs_pairs <- lapply(na.omit(un_tr_pairs[,-(1:5)]), function(i) wilcox.test(i ~ un_tr_pairs$Treatment, paired=TRUE)$p.value )
diffs_pairs <- unlist(diffs_pairs)
diffs_pairs1 <- diffs_pairs[diffs_pairs < 0.05]


# Figure 5B
# Volcano plot
unt <- un_tr_pairs[un_tr_pairs$Treatment == "Pre",c(6:138)]
unt1 <- as.data.frame(sapply(unt, mean))
colnames(unt1)[1] <- "Untreated"
treat <- un_tr_pairs[un_tr_pairs$Treatment == "Post",c(6:138)]
treat1 <- as.data.frame(sapply(treat, mean))
colnames(treat1)[1] <- "Treated"
df <- merge(unt1, treat1, by="row.names")
df$Fold <- df$Treated/df$Untreated
df$Log2FC <- log2(df$Fold)
pvals <- as.data.frame(unlist(diffs_pairs))
colnames(pvals)[1] <- "pvalue"
df2 <- merge(df, pvals, by.x="Row.names", by.y="row.names")
df2$sig <- ifelse(df2$pvalue < 0.05, df2$sig <- "p<0.05", df2$sig <- "ns")
ggplot(df2, aes(x = Log2FC, y = -log10(pvalue))) +
  geom_point(aes(color = sig)) +
  scale_color_manual(values = c("grey", "red")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  geom_text_repel(
    data = subset(df2, pvalue <= 0.05), 
    aes(label = Row.names),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    segment.colour = "transparent") +
  labs(y="-Log10(pValue)", x="Log2(Treated/ untreated)")


# Figure 5C
# Paired dotplots
melt_olink <- reshape2::melt(olinkX_ID, id.vars = c("Sample", "Group2", "Batch2", "Treatment","Subject"))
melt_olink <- subset(melt_olink, variable %in% names(diffs_pairs1))
melt_olink <- subset(melt_olink, Subject %in% MISCpairs)
melt_olink$Treatment <- factor(melt_olink$Treatment, levels = c("Pre", "Post"))
melt_olink$variable <- factor(melt_olink$variable, levels = c("CXCL10", "EN-RAGE", "HEXIM1", "PSIP1", "CCL25", "ITGA11", "TNFB"))
ggplot(melt_olink, aes(x=Treatment, y=value)) + theme_minimal() + 
  geom_jitter(shape=16, position=position_jitter(0.02), size=3, aes(colour=Treatment)) +
  geom_line(aes(group=Subject), position=position_dodge(width=0.1)) +
  facet_wrap(~ variable, scales = "free", ncol = 4) + 
  scale_color_manual(values=c("#A5AA99", "lightgrey")) +
  theme(axis.text.x = element_text(size = 10), legend.title = element_blank()) + labs(x=NULL, y="NPX") 




sessionInfo()
#R version 3.6.2 (2019-12-12)
#Platform: x86_64-apple-darwin15.6.0 (64-bit)
#Running under: macOS Catalina 10.15.6

#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] limma_3.40.6       reshape2_1.4.4     OlinkAnalyze_1.0.2 readxl_1.3.1       forcats_0.5.0      dplyr_1.0.0        purrr_0.3.4       
#[8] readr_1.3.1        tidyr_1.1.0        tibble_3.0.2       tidyverse_1.3.0    stringr_1.4.0      openxlsx_4.1.5     lmerTest_3.1-2    
#[15] lme4_1.1-23        Matrix_1.2-18      ggrepel_0.8.2      ggfortify_0.4.10   ggplot2_3.3.2      emmeans_1.4.8      car_3.0-8         
#[22] carData_3.0-4      broom_0.7.0       



