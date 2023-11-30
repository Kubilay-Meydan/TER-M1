rule salmon_quant_pe:
    input:
        fast1 = lambda wildcards: get_processed_fastq(wildcards.sample, pair=1),
        fast2 = lambda wildcards: get_processed_fastq(wildcards.sample, pair=2),
        index = os.path.join(TXOME_DIR, "seq.bin")

    output:
        os.path.join(OUTPUT_DIR, "{sample}", "quant.sf")

    params:
        index_dir = TXOME_DIR,
        output_dir = os.path.join(OUTPUT_DIR, "{sample}"),
        libtype = "A",
        gtf = get_gtf(SPECIES),
        extra_params = return_parsed_extra_params(config["extra_salmon_parameters"]),
        threads = cluster_dict["salmon_quant_pe"]["smp"]

    # threads: 4
    conda:
        "../env/align.yaml"

    shell:
        """
        salmon quant \
        --index {params.index_dir} \
        --libType {params.libtype} \
        --mates1 {input.fast1} \
        --mates2 {input.fast2} \
        --geneMap {params.gtf} \
        --threads {params.threads} \
        {params.extra_params} \
        -o {params.output_dir} \
        """
