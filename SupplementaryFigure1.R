#Created by Giuseppe Rubens Pascucci (pascucci.1479790@studenti.uniroma1.it)
####################################################################
# pip install mofapy2 (on terminal)
# devtools::install_github("bioFAM/MOFA2/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"))
# https://github.com/bioFAM/MOFA2/blob/master/MOFA2/vignettes/getting_started.md
# https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/downstream_analysis.html
# (Optional) set up reticulate connection with Python
# reticulate::use_python("/Users/ricard/anaconda3/envs/base_new/bin/python", required = T)

library(MOFA2)
library(ggplot2)
library(ggbeeswarm)
library(GGally)
library(tidyr)
library(ggpubr)

mycols <- c(Kawasaki = '#bd377f', 'CoV2+Kawasaki' = '#9da290', 'CoV2+' = '#d7802f', Controls = '#3780b2')


# ================================================== DATABASE LOAD ==================================================

data <- read.table("DataIN/Dataset_FACS+Prot.txt", sep = "\t", header = T, check.names=F, stringsAsFactors = F)
data <- as.matrix(data)
rownames(data) <- data[,1]
Group <- data[,2]
Age <- data[,3]
Sex <- data[,4]
data <- data[,-c(1:4)]
class(data) <- "numeric"

prot <- data[, -c(1:19)]
freq <- data[, -c(20:140)]


# ================================================== DATA LOAD ==================================================

ldata <- list(FACS = t(freq), Proteomics = t(prot))


# ================================================== CREATE MOFA OBJECT ==================================================

MOFAobject <- create_mofa(ldata)

# Visualise data structure
plot_data_overview(MOFAobject)


# ================================================== DEFINE OPTIONS ==================================================

# Data options
data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views <- TRUE

# Model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 10

# Training options
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "medium"
train_opts$seed <- 42

# ================================================== PREPARE MOFA OBJECT ==================================================

MOFAobject <- prepare_mofa(MOFAobject, data_options = data_opts, model_options = model_opts, training_options = train_opts) # stochastic_options = stochastic_opts


# ================================================== TRAIN THE MODEL ==================================================

# R requires the package reticulate to communicate with Python. 
# https://rstudio.github.io/reticulate/
library(reticulate)

# Using a conda enviroment called myR
use_condaenv("myR", required = TRUE)

outfile <- paste0(getwd(),"/MOFA/model.hdf5")
MOFAmodel <- run_mofa(MOFAobject, outfile)



# ================================================== LOAD TRAINED MODEL AND ADD METADATA ==================================================

model <- load_model(paste0(getwd(),"/MOFA/model.hdf5"))

sample_metadata <- model@samples_metadata
head(sample_metadata, n=3)
sample_metadata$condition <- Group
sample_metadata$age <- Age
sample_metadata$sex <- Sex
sample_metadata

samples_metadata(model) <- sample_metadata
head(model@samples_metadata, n=3)


# ================================================== VARIANCE DECOMPOSITION ==================================================

# Total variance explained per view and group
head(model@cache$variance_explained$r2_total[[1]]) # group 1

# Variance explained for every factor in per view and group
head(model@cache$variance_explained$r2_per_factor[[1]]) # group 1

plot_variance_explained(model, x="view", y="group") +
  theme(axis.text.x = element_text(color="black", angle=50, vjust=1.4, hjust= 1.2, size = 12), axis.text.y = element_text(color="white"))
ggsave(("MOFA/Variance Decomposition v1.pdf"), width = 8, height = 5)

plot_variance_explained(model, x="view", y="factor") +
  theme(axis.text.x = element_text(angle=50, vjust=1, hjust= 1, size = 15), axis.text.y = element_text(size = 15), legend.title = element_text(size = 13))
ggsave(("MOFA/Variance Decomposition v2.pdf"), width = 4, height = 9)

plot_variance_explained(model, x="view", y="factor", plot_total = T)[[2]] +
  theme(axis.text.x = element_text(angle=50, vjust=1, hjust= 1, size = 19), axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 15))
ggsave(("MOFA/Variance Decomposition Tot.pdf"), width = 3, height = 7)



# ================================================== VISUALIZATION OF SINGLE FACTORS ==================================================

# Get factor values
Z <- get_factors(model, factors = 1:10 , groups = "all", as.data.frame=TRUE)
Z$factor <- as.factor(Z$factor)
df <- merge(Z,  cbind(sample = rownames(data), condition = Group, sex = Sex, age = Age), by="sample")
df$age <- as.numeric(as.character(df$age))

# Beeswarm Plots v1
ggplot(df, aes(x=factor, y=value, col=condition, shape=sex)) + facet_wrap(~factor, nrow=1, scales="free_x") +
  labs(x="", y="Factor value") + geom_quasirandom(size = 5, dodge.width = 0, alpha = 0.8, method = "quasirandom", width = 0.4) + 
  scale_color_manual(values=mycols) + theme_classic() + geom_hline(yintercept=0, linetype="dashed", size=0.4, alpha=0.5) +
  scale_shape_manual(values=c(20, 18)) +
  theme(text = element_text(size=18), panel.border = element_rect(color="black", size=0.2, fill=NA), 
        strip.background = element_rect(colour = "black", size=0.3), panel.spacing = unit(0,"lines"),
        axis.line = element_line(colour = 'grey', size = 0.2), axis.text.x = element_text(color="white"),
        axis.ticks.x = element_line(color = "white"))
ggsave(("MOFA/Beeswarm plots v1.pdf"), width = 14, height = 6)

# Beeswarm Plots v2
plot_factor(model, factor = 1:10, color_by = "condition", shape_by = "sex", dot_size = 4, dot_alpha = 0.8, scale = T, dodge = F) + 
  scale_fill_manual(values=mycols) +
  theme(text = element_text(size=18))
ggsave(("MOFA/Beeswarm plots v2.pdf"), width = 14, height = 6)

# Violin Plots v1
ggplot(df, aes(x=factor, y=value, col=condition)) + facet_wrap(~factor, nrow=1, scales="free_x") +
  labs(x="", y="Factor value") + geom_quasirandom(size = 1.5, dodge.width = 3, alpha = 0.8, method = "quasirandom", width = 0) + 
  geom_violin(alpha= 0.5, trim = TRUE, scale = "width", position = position_dodge(width=3), show.legend = FALSE) +
  scale_color_manual(values=mycols) + theme_classic() + geom_hline(yintercept=0, linetype="dashed", size=0.4, alpha=0.5) +
  scale_shape_manual(values=c(20, 18)) +
  theme(text = element_text(size=18), panel.border = element_rect(color="black", size=0.2, fill=NA), 
        strip.background = element_rect(colour = "black", size=0.3), panel.spacing = unit(0,"lines"),
        axis.line = element_line(colour = 'grey', size = 0.2), axis.text.x = element_text(color="white"),
        axis.ticks.x = element_line(color = "white"))
ggsave(("MOFA/Violin plots v1.pdf"), width = 14, height = 6)

# Violin Plots v2
plot_factor(model, factor = 1:10, color_by = "condition", dot_size = 1.5, dodge = T, legend = T, add_violin = T, violin_alpha = 0.3) + 
  scale_fill_manual(values=mycols) + theme(text = element_text(size=18)) 
ggsave(("MOFA/Violin plots v2.pdf"), width = 14, height = 6)



# ================================================== VISUALIZATION OF COMBINATION OF FACTORS ==================================================

# spread over factors
df2 <- spread(df, key="factor", value="value")

ggpairs(df2, columns = c(6:10),
        lower = list(continuous=GGally::wrap("points", size=2, alpha = 0.6)), 
        diag = list(continuous=GGally::wrap("densityDiag", alpha = 0.5, size = 0.2)), 
        upper = list(continuous=GGally::wrap("points", size=2, alpha = 0.6)),
        mapping = ggplot2::aes(color=condition)) + 
  scale_fill_manual(values=mycols) + scale_color_manual(values=mycols) + 
  theme(text = element_text(size=16), panel.grid.major = element_blank(),
        panel.border = element_rect(color="black", size=0.3, fill=NA),
        panel.background = element_rect(fill = "white"), strip.background = element_rect(color = "black", size=0.3, fill = "gray94"))
ggsave(("MOFA/Combinations of factors (1-5).pdf"), width = 7, height = 7)



# ================================================== VISUALIZATION OF FEATURE WEIGHTS ==================================================

plot_top_weights(object = model, view = "FACS", factor = 5, nfeatures = 10) + 
  theme(text = element_text(size=15), strip.background = element_rect(color = "black", fill = "gray94"))
ggsave(("MOFA/Feature weights F5 (FACS).pdf"), width = 10, height = 4)
  
plot_top_weights(object = model, view = "Proteomics", factor = 5, nfeatures = 10) + 
  theme(text = element_text(size=15), strip.background = element_rect(color = "black", fill = "gray94"))
ggsave(("MOFA/Feature weights F5 (Proteomics).pdf"), width = 10, height = 4)


plot_top_weights(object = model, view = "FACS", factor = 1, nfeatures = 10) + 
  theme(text = element_text(size=15), strip.background = element_rect(color = "black", fill = "gray94"))
ggsave(("MOFA/Feature weights F1 (FACS).pdf"), width = 10, height = 4)

plot_top_weights(object = model, view = "Proteomics", factor = 1, nfeatures = 10) + 
  theme(text = element_text(size = 15), strip.background = element_rect(color = "black", fill = "gray94"))
ggsave(("MOFA/Feature weights F1 (Proteomics).pdf"), width = 7, height = 4)


