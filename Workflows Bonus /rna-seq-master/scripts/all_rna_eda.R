#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
txi_rdata = args[1]
pca_plot_pdf = args[2]
out_rdata = args[3]
factor_str = args[4]
library_tsv = args[5]

## txi_rdata = "test/analysis/all_txi.rdata"
## library_tsv = "test/inputs/libraries.tsv"
## factor_str = "run group"
## out_rdata = "test/analysis/eda.rdata"
## pca_plot_pdf = "test/results/all_pca.pdf"

library(cowplot)

library(ggrepel)
library(tidyverse)

load(txi_rdata)

library(edgeR)

counts = txi$counts

# Obtaining per-observation scaling factors for length, adjusted to avoid
# (see https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#edgeR)
norm_mat = txi$length
norm_mat = norm_mat/exp(rowMeans(log(norm_mat)))
norm_counts = counts/norm_mat

# Get effective library sizes from scaled counts
eff_lib = calcNormFactors(norm_counts) * colSums(norm_counts)
norm_mat = sweep(norm_mat, 2, eff_lib, "*")
norm_mat = log(norm_mat)

# Creating a DGEList object for use in edgeR.
y = DGEList(counts)
y = scaleOffset(y, norm_mat)

libraries = read_tsv(library_tsv)
factor_vec = unlist(strsplit(factor_str, " "))

formula = as.formula(paste("~ ", paste(factor_vec, collapse = "+")))
formula
design = model.matrix(formula, libraries)
design

# filtering using the design information
## design <- model.matrix(~condition, data = sampleTable)
keep <- filterByExpr(y, design)
y <- y[keep, ]

logCPM <- cpm(y, prior.count=2, log=TRUE, offset = y$offset)

pca = prcomp(t(logCPM))

make_pca_plots = function(in_pca, full_libs){
  pve_pc1=round(100*summary(in_pca)$importance[2,1])
  pve_pc2=round(100*summary(in_pca)$importance[2,2])
  pca_plot = as.data.frame(in_pca$x) %>%
    rownames_to_column(var = "library") %>%
    left_join(libraries, by = "library") %>%
    ggplot(., aes(x = PC1, y = PC2, color = get(tail(factor_vec, n= 1)), label = library)) +
    geom_point(size = 4) +
    geom_text_repel() +
    xlab(paste("PC1, ", pve_pc1, "% variance explained", sep ="")) +
    ylab(paste("PC2, ", pve_pc2, "% variance explained", sep ="")) +
    scale_color_discrete(name = paste0(tail(factor_vec, n=1))) +
    theme_cowplot() +
    theme(legend.position = "bottom")
  return(pca_plot)
}

pca_plot = make_pca_plots(pca, libraries)
# Note this resembles plotMDS(y, gene.selection = "common")

save(design, formula, logCPM, pca, pca_plot, y, file = out_rdata)

save_plot(pca_plot, file = pca_plot_pdf)
