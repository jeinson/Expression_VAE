# This script will manipulate the covariates from the VAE model to perfectly
# match the covariates input for qtltools

setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/random/ML_genomics")
rm(list = ls())
library(tidyverse)

peer_factors <- read_tsv("Muscle_Skeletal.v8.covariates.txt")

vae_encodings <- read.csv("VAE_sampled_encodings_mod8.tsv", 
                          sep = "\t", row.names = 1)
vae_encodings <- t(vae_encodings)
rownames(vae_encodings) <- paste0("Encoding", 1:60)
vae_encodings <- as.data.frame(vae_encodings)
vae_encodings <- rownames_to_column(vae_encodings, "ID")

# Compile the output matrix
# 
# Start with taking the genetic PCs
out <- peer_factors[1:5,]
out <- rbind(out, vae_encodings[,c("ID", colnames(out)[-1])])
out <- rbind(out, peer_factors[66:68,]) # Add the sex, pcr, and sequencing covariates

write_tsv(out, "Muscle_Skeletal_VAE_covariates.txt")

# Generate some totally random covariates to map eQTLs with
hist(unlist(vae_encodings[,-1]))
muhat = mean(unlist(vae_encodings[,-1]))
sdhat = sd(unlist(vae_encodings[,-1]))

n_tot = length(unlist(vae_encodings[,-1]))

control_encodings <- rnorm(n_tot, muhat, sdhat)
control_encodings <- matrix(control_encodings, nrow = nrow(vae_encodings))
control_encodings <- data.frame(vae_encodings$ID, control_encodings)
colnames(control_encodings) <- colnames(vae_encodings)

control_out <- peer_factors[1:5,]
control_out <- rbind(control_out, control_encodings[,c("ID", colnames(control_out)[-1])])
control_out <- rbind(control_out, peer_factors[66:68,]) # Add the sex, pcr, and sequencing covariates

write_tsv(control_out, "FAKE_Muscle_Skeletal_encodings.txt")
