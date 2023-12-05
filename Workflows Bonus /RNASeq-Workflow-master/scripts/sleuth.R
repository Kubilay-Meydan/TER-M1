library("sleuth")

# Specify where kallisto results are stored.
sample_ids = t(read.delim(file.path(snakemake@params["samples_info"]), header = TRUE, sep = "\t"),["sample_name"])

# Load table that describes samples and conditions, and links those to kallisto output.
kal_dirs <- file.path(snakemake@params["kal_dirs"], sample_name)
s2c <- read.table(file.path(snakemake@params["samples_info"]), header=TRUE, sep="\t", stringsAsFactors=FALSE)
s2c <- dplyr::mutate(s2c, path=kal_dirs)

# Construct sleuth object.
so <- sleuth_prep(s2c, extra_bootstrap_summary=TRUE)

# Fit full model.
so <- sleuth_fit(so, ~condition, "full")

# Fit reduced model.
so <- sleuth_fit(so, ~1, "reduced")

# Apply likelihood-ratio test and make sleuth table.
so <- sleuth_lrt(so, "reduced", "full")
sleuth_table <- sleuth_results(so, "reduced:full", "lrt", show_all=FALSE)

# Get the significant genes (p-value <= 0.05).
sleuth_significant <- dplyr::filter(sleuth_table, p_val <= 0.05)
write(t(sleuth_significant), file = file.path(snakemake@params["significant_genes"]), sep="\t")

# Generate and save plots.
# PCA
pca_plot <- plot_pca(so, color_by="condition")
ggplot2::ggsave(filename=snakemake@params["pca_plot"], plot=pca_plot)

# MA
ma_plot <- plot_ma(so)
ggplot2::ggsave(filename=snakemake@params["ma_plot"], plot=ma_plot)

# Volcano
volcano_plot <- plot_volcano(so)
ggplot2::ggsave(filename=snakemake@params["volcano_plot"], plot=volcano_plot)

# Sample heatmap
sample_heatmap <- plot_sample_heatmap(so)
ggplot2::ggsave(filename=snakemake@params["sample_heatmap"], plot=sample_heatmap)
