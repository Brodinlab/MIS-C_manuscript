######################################
#                PCA                 #
######################################
library(ggfortify)
library(devtools)

KD <- read.csv('Olink_KD_CoV2_NA_30cutoff_minTreated.csv', sep = ';')
df <- KD[,-c(1:2)]
head(df)

pca_res <- prcomp(df, scale. = TRUE)

autoplot(pca_res, data = KD, colour = 'KD_COV') 
pairs(pca_res$loadings)

library(factoextra)
res.var <- get_pca_var(pca_res)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 

fviz_pca_var(pca_res, select.var = list(contrib = 20), repel=TRUE, 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")) 

fviz_pca_var(pca_res, alpha.var = "contrib")

fviz_pca_ind(pca_res)

# Contributions of variables to PC1
fviz_contrib(pca_res, choice = "var", axes = 1, top = 20)
# Contributions of variables to PC2
fviz_contrib(pca_res, choice = "var", axes = 2, top = 20)
