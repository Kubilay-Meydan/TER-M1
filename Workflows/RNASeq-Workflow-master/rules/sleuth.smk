rule sleuth:
    input:
        expand(config["kallisto_out"] + "/{dataset}/abundance.h5", dataset = list(samples_info.index)),
        expand(config["kallisto_out"] + "/{dataset}/abundance.tsv", dataset = list(samples_info.index)),
        expand(config["kallisto_out"] + "/{dataset}/run_info.json", dataset = list(samples_info.index)),
        config["samples_info"]
    output:
        config["sleuth_out"] + "sleuth_output.tsv"
    conda:
        "../envs/sleuth.yaml"
    params:
        sample_tsv = config["samples_info"]
        kal_dirs = config["kallisto_out"]
        significant_genes = congif["sleuth_out"]
        pca_plot = config["sleuth_out"] + "/pca_plot.pdf"
        ma_plot = config["sleuth_out"] + "/ma_plot.pdf"
        volcano_plot = config["sleuth_out"] + "/volcano_plot.pdf"
    script:
        "../scripts/sleuth.R"
