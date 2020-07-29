# Plasma protein analyses adult COVID19 vs MIS-C and Kawasaki cohorts (Figure 2)

######################################
#           Normalization            #
######################################
# install.packages("devtools")
devtools::install_github(repo ='Olink-Proteomics/OlinkRPackage/OlinkAnalyze')

library(OlinkAnalyze)
# Load olink data
May <-read_NPX(filename= 'Petter_Brodin_Covid-19_Plate2_NPX_LOD.xlsx')
July <- read_NPX(filename = "OlinkJuly1.xlsx")
July$NPX[July$NPX <= July$LOD] <- July$LOD[July$NPX <= July$LOD] # substitute NPX values under LOD for LOD
# correct values for samples that were diluted
July$unnormalized <- 2^(July$NPX)
dilutedby2 <- c("CACTUS 004", "CACTUS 031", "CACTUS 032", "CACTUS 018", "CACTUS 023", "CACTUS 024", "CACTUS 017")
July$dilutionfactor <- ifelse(July$SampleID %in% dilutedby2,
                              July$dilutionfactor <- 2,     
                              July$dilutionfactor <- 1 )
July$unnormalized2 <- (July$unnormalized)*(July$dilutionfactor)
July$NPX2 <- log2(July$unnormalized2)
#new df for july
newJuly <- July
newJuly$NPX <- newJuly$NPX2  
newJuly <- newJuly[,-c(14:17)]
rm(July, dilutedby2)
# identify bridge samples
bridge_samples <- intersect(x = newJuly$SampleID, y = May$SampleID)
# bridge normalization (using bridge samples common in both runs)
bridge_normalized_data <- olink_normalization(df1 = May, df2 = newJuly, overlapping_samples_df1 = bridge_samples)

# cast df, this will average out proteins that were run between two panels
olink <- bridge_normalized_data[,c(1,5,12)]
library(reshape2)
olink <- dcast(olink, SampleID ~ factor(Assay, levels = unique(olink$Assay)), fun.aggregate=function(i) mean(i, na.rm=TRUE))
rownames(olink) <- olink$SampleID
olink <- olink[,-1]

write.xlsx(olink, 'Olink_Kaw_MISC.xlsx', row.names=TRUE)

######################################
#             Z-score                #
######################################
# Data tables re-formatted as columns (proteins) and rows (samples)
RUN12 <- read.csv('Olink_all_run12.csv', sep = ';')
df <- RUN12[,-1]
RUN12_Z <- apply(df, 1, scale)
write.csv(t(RUN12_Z), 'Olink_all_run12Z.csv')

# File 'Petter_Brodin_Covid-19_Plate1_NPX_LOD.xlsx' where Adult CoV2+ cases located, filtered by 30% cutoff and 
# matching proteins between this plate and the normalized set of children cases
RUN3 <- read.csv('Olink_all_run3.csv', sep = ';')
df <- RUN3[,-1]
RUN3A_Z <- apply(df, 1, scale)
write.csv(t(RUN3A_Z), 'Olink_all_run3AZ.csv')

a <- read.csv('Olink_all_run12Z.csv', sep=';')
b <- read.csv('Olink_all_run3AZ.csv', sep=';')
KD <- rbind(a,b)

######################################
#                PCA                 #
######################################
library(ggfortify)
library(devtools)

df <- KD[,-c(1:2)]
head(df)

pca_res <- prcomp(df, scale. = TRUE)

#PCA plot PC1 and PC2
autoplot(pca_res, data = KD, colour = 'KD_COV') 
pairs(pca_res$loadings)

library(factoextra)
res.var <- get_pca_var(pca_res)

#PC contribution plot
fviz_pca_var(pca_res, select.var = list(contrib = 20), repel=TRUE, 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")) 

fviz_pca_var(pca_res, alpha.var = "contrib")

fviz_pca_ind(pca_res) 

# Contributions of variables to PC1
fviz_contrib(pca_res, choice = "var", axes = 1, top = 20)
# Contributions of variables to PC2
fviz_contrib(pca_res, choice = "var", axes = 2, top = 20)
