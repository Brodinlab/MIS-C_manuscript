# The Immunology of Multisystem Inflammatory Syndrome in Children with COVID-19
# Figure 6

# Load preprocessed data
virscan <- read.csv("https://ki.box.com/shared/static/sze052yaspnhe2czp8j6qgniix6tc3xh.csv", row.names = 1, stringsAsFactors = FALSE)
label <- read.csv("https://ki.box.com/shared/static/9c84p9ue2gffvys2vttbdrphvdfm59ch.csv", row.names = 1, stringsAsFactors = FALSE)
m_Vir <- merge(label, virscan, by="row.names")
m_Vir$Label <- factor(m_Vir$Label, levels = c("Healthy","MIS-C", "COV2", "Kawasaki"))
colnames(m_Vir)[1] <- "ID"

# Figure 6C
library(stringr)
corona <- c(colnames(m_Vir)[str_detect(colnames(m_Vir), "corona")], colnames(m_Vir)[str_detect(colnames(m_Vir), "cov")])
corona <- m_Vir[,c("ID", "Label", corona)]
library(reshape2)
corona <- reshape2::melt(corona, id.vars=c("ID", "Label"))
corona$Label <- factor(corona$Label, levels = c("Healthy", "COV2", "MIS-C", "Kawasaki"))
library(ggplot2)
ggplot(corona, aes(y=value, x=variable)) +
  geom_jitter(shape=16, height = 0.1, aes(color=Label, size=value)) +
  theme_minimal() + 
  theme(legend.title = element_blank()) +
  scale_color_manual(values=c("#2F8AC4", "#E48725", "#A5AA99", "#CD3A8E", "black")) +
  facet_wrap(~ Label, ncol = 1, scales = "free") + labs(x=NULL, y="VirScore") +
  ylim(c(-1,7)) + coord_flip()



sessionInfo()
#R version 3.6.2 (2019-12-12)
#Platform: x86_64-apple-darwin15.6.0 (64-bit)
#Running under: macOS Catalina 10.15.6

#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] ggplot2_3.3.2  reshape2_1.4.4 stringr_1.4.0 
