# Analysis of eQTLs fit with VAE covariate parameters
# 
# ML Genomics Final Project
# Jonah Einson

rm(list = ls())
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/random/ML_genomics")
library(tidyverse)

# Read in eqtl results from the vae method
eqtl_result_paths <- list.files("12.07.19_eqtl_results", full.names = T)

eqtl_result_paths <- 
  eqtl_result_paths[order(nchar(eqtl_result_paths), eqtl_result_paths)]

vae_eqtls <- data.frame() 
for(i in eqtl_result_paths){
  vae_eqtls <- 
    rbind(
      vae_eqtls, 
      read_delim(i, delim = " " ,col_names = F)
    )
}
qtl_colnames  <- 
  c("gene_id",     #1
    "gene_chr",    #2
    "gene_start",  #3
    "gene_end",    #4 
    "strand",      #5 
    "nvar",        #6
    "tss_distance",#7 
    "variant_id",  #8
    "chr",         #9
    "pos",         #10  
    "end",         #11
    "true_df",     #12
    "dummy",       #13
    "beta_shape1", #14
    "beta_shape2", #15
    "pval_nominal",#16  
    "slope",       #17 
    "pval_true_df",#18
    "pval_beta"    #19)
  )
colnames(vae_eqtls) <- qtl_colnames

# Adjust for multiple comparisons with the qvalue package
library(qvalue)
qvals <- qvalue(vae_eqtls$pval_beta)
hist(qvals)
vae_eqtls$qval <- qvalue(vae_eqtls$pval_beta)$qval 

# Read in the "official" eQTL results
eqtl_path <- "/gpfs/commons/datasets/controlled/GTEx/dbgap_restricted/data/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL/Muscle_Skeletal.v8.egenes.txt.gz"
gtex_eqtls <- read_tsv(eqtl_path)

official_qvals <- qvalue(gtex_eqtls$pval_beta)
hist(official_qvals)

# Read in eQTLs mapped without covariates
eqtl_result_paths <- list.files("12.08.19_eqtl_results", full.names = T)

eqtl_result_paths <- 
  eqtl_result_paths[order(nchar(eqtl_result_paths), eqtl_result_paths)]

noco_eqtls <- data.frame() 
for(i in eqtl_result_paths){
  noco_eqtls <- 
    rbind(
      noco_eqtls, 
      read_delim(i, delim = " " ,col_names = F)
    )
}
colnames(noco_eqtls) <- qtl_colnames

noco_qvals <- qvalue(noco_eqtls$pval_beta)
hist(noco_qvals)
noco_eqtls$qval <- noco_qvals$qvalues

# Read in eQTLs mapped with dummy covariates
eqtl_result_paths <- list.files("12.08.19_fake_eqtl_results/", full.names = T)
eqtl_result_paths <- 
  eqtl_result_paths[order(nchar(eqtl_result_paths), eqtl_result_paths)]

fake_eqtls <- data.frame() 
for(i in eqtl_result_paths){
  fake_eqtls <- 
    rbind(
      fake_eqtls, 
      read_delim(i, delim = " " ,col_names = F)
    )
}
colnames(fake_eqtls) <- qtl_colnames

fake_qvals <- qvalue(fake_eqtls$pval_beta)
hist(fake_qvals)
fake_eqtls$qval <- fake_qvals$qvalues


# Let's do some comparison
dim(gtex_eqtls)
dim(vae_eqtls)

# Make sure we're working with the same set of genes
shared_genes <- intersect(gtex_eqtls$gene_id, vae_eqtls$gene_id)

sum(gtex_eqtls$qval < .05)
sum(vae_eqtls$qval < .05)
gtex_eqtls <- gtex_eqtls[match(shared_genes, gtex_eqtls$gene_id),]
vae_eqtls <- vae_eqtls[match(shared_genes, vae_eqtls$gene_id),]
noco_eqtls <- noco_eqtls[match(shared_genes, noco_eqtls$gene_id),]
fake_eqtls <- fake_eqtls[match(shared_genes, fake_eqtls$gene_id),]


npval_cor <- cor(gtex_eqtls$pval_beta, vae_eqtls$pval_beta, method = "spearman")

# Plot the beta p-values against each other
plt_tbl = data.frame(gtex_beta_pvals = gtex_eqtls$qval, 
                     vae_beta_pvals = vae_eqtls$qval)

beta_peer_pvalue_plot <- 
  ggplot(plt_tbl, aes(gtex_beta_pvals, vae_beta_pvals)) + 
  geom_point(alpha = .1, color = rgb(21, 118, 178, maxColorValue = 255)) +
  geom_rug(alpha = .005) +
  theme_classic() +
  ylab("VAE covariate eQTL q-values") + 
  xlab("Gold Standard eQTL q-values")

png("final_figs/peer_qvalue_plot.png", width = 5, height = 4, units = "in", res = 400)
beta_peer_pvalue_plot
dev.off()

 # What is the overlap of significant p-values?
vae_egenes <- vae_eqtls$gene_id[vae_eqtls$qval < .05]
gtex_egenes <- gtex_eqtls$gene_id[gtex_eqtls$qval < .05]
noco_egenes <- noco_eqtls$gene_id[noco_eqtls$qval < .05]
fake_egenes <- fake_eqtls$gene_id[fake_eqtls$qval < .05]

require(gplots)
gene_list = list(PEER = gtex_egenes, VAE = vae_egenes, None = noco_egenes, Fake = fake_egenes)
venn(list(PEER = gtex_egenes, VAE = vae_egenes, None = noco_egenes, Fake = fake_egenes))

library(UpSetR)
overlap_plt <- 
  upset(fromList(gene_list), 
        order.by = "freq", 
        sets.x.label = "Total eGenes", mainbar.y.label = "n intersecting eGenes",
        matrix.color = rgb(21, 118, 178, maxColorValue = 255))

png("final_figs/overlap-venn-upset-plt.png", width = 5, height = 3, units = "in", res=400)
overlap_plt
dev.off()

library(GeneOverlap)
go.obj <- newGeneOverlap(gtex_egenes, 
                         vae_egenes, 
                         genome.size = length(shared_genes) # The number of co-covered genes
)

go.obj <- testGeneOverlap(go.obj)
go.obj # Clearly there's some overlap here. Not too sure what to make of it though



