
rule salmon_quant_se:
    input:
        fast1 = lambda wildcards: get_processed_fastq(wildcards.sample, pair=1),
        index = os.path.join(TXOME_DIR, "seq.bin")

    output:
        os.path.join(OUTPUT_DIR, "{sample}", "quant.sf")

    params:
        index_dir = TXOME_DIR,
        output_dir = os.path.join(OUTPUT_DIR, "{sample}"),
        libtype = "A",
        gtf = get_gtf(SPECIES),
        extra_params = return_parsed_extra_params(config["extra_salmon_parameters"]),
        threads = cluster_dict["salmon_quant_se"]["smp"]

    conda:
        "../env/align.yaml"

    shell:
        """
        salmon quant \
        --index {params.index_dir} \
        --libType {params.libtype} \
        -r {input.fast1} \
        --geneMap {params.gtf} \
        --threads {params.threads} \
        {params.extra_params} \
        -o {params.output_dir} \
        """
