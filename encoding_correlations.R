# Look for correlations between the VAE encodings and known features in the
# dataset. 
# 
# Jonah Einson
# Genomics ML final Project

rm(list = ls())
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/random/ML_genomics")
library(tidyverse)

encodings <- read.csv("VAE_sampled_encodings_mod8.tsv", 
                      sep = "\t", row.names = 1,
                      col.names = c(paste0("Encoding", 1:60))
                      )

metadata <- read.csv("/gpfs/commons/datasets/controlled/GTEx/dbgap_restricted/data/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt", 
                     sep = "\t",
                     row.names = 1)

# Does this look like anything off the bat?
plot(encodings$Encoding1, encodings$Encoding2)

# First let's see if any of these encoding correlates with sex
sex <- metadata[rownames(encodings),]$SEX
sex_corrs <- 
  apply(encodings, 2, function(x){
    cor(x, sex, method = "spearman")
  })

plot(sex_corrs, type = "l")

x1 <- names(sort(sex_corrs)[60])
x2 <- names(sort(sex_corrs)[59])

plot(encodings[[x1]], encodings[[x2]], 
     xlab = x1, ylab = x2, pch = 19,
     col = sex)


library(corrplot)
encoding_autocorrelation <- cor(encodings[,-1])
corrplot(cor(encodings[,-1]), method = "color", tl.cex = .5, tl.col = "black")

# Plot a couple of the encodings which are strongly correlated
corr_features <- 
  which(
    abs(encoding_autocorrelation) > .6 & abs(encoding_autocorrelation) < 1, 
    arr.ind = T
  )

plot(encodings[, corr_features[7,]])

############ See if the VAE encodings correlate with PEER factors ##############
peer_covars <- read_tsv("Muscle_Skeletal.v8.covariates.txt")
peer_covars <- peer_covars[startsWith(peer_covars$ID, "Infer"),]
peer_covars %<>% as.data.frame %>% column_to_rownames("ID") 

vae_correlation_matrix <- matrix(nrow = 60, ncol = 60, 
                                 dimnames = list(rownames(peer_covars), 
                                                 colnames(encodings)[1:60]))

for(i in 1:60){
  for(j in 1:60){
    vae_correlation_matrix[i,j] <- cor(unlist(peer_covars[i,]), encodings[,j], method = "spearman")
  }
}

sum(abs(vae_correlation_matrix))

vae_correlation_matrix %>%
  as.data.frame %>%
  rownames_to_column("PEER") %>%
  gather(vae_features, corr, - PEER) %>%
  mutate(PEER = factor(PEER, levels = paste0("InferredCov", 1:60))) %>%
  mutate(vae_features =  factor(vae_features, levels = rev(paste0("Encoding", 1:60)))) ->
  vae_correlation_matrix_tidy

peer_vae_correlations_plt <- 
  ggplot(vae_correlation_matrix_tidy, aes(PEER, vae_features, fill = corr)) +
  geom_tile() + 
  theme_minimal() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.position = c(.8, .8) ) +
  ylab("VAE Encoding")

png("final_figs/vae_peer_correlaton_final.png", width = 7, height = 7, units = "in", res = 400)
peer_vae_correlations_plt
dev.off()

barplot(colSums(abs(t(vae_correlation_matrix))), 
        col = rgb(21, 118, 178, maxColorValue = 255), space = 0, border = NA, 
        xaxt = "n")

barplot(rowSums(abs(t(vae_correlation_matrix))), 
        col = rgb(21, 118, 178, maxColorValue = 255), space = 0, border = NA, 
        xaxt = "n")
