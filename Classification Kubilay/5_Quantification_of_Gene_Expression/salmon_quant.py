rule salmon_quant:
    input:
        fast1 = lambda wildcards: get_processed_fastq(wildcards.sample, pair=1),
        fast2 = lambda wildcards: get_processed_fastq(wildcards.sample, pair=1),
        index = os.path.join(scallop_outdir, "extended_transcriptome/seq.bin"),
        scallop_ref = os.path.join(scallop_outdir,"scallop.tx_gene.tsv")
    output:
        os.path.join(scallop_outdir, "{sample}", "quant.sf")
    params:
        salmon = config["salmon_path"],
        index_dir = os.path.join(scallop_outdir, "extended_transcriptome/"),
        output_dir = os.path.join(scallop_outdir, "{sample}"),
        libtype = get_salmon_strand(config["feature_counts_strand_info"]),
        extra_params = return_parsed_extra_params(config["extra_salmon_parameters"])
    threads: 4
    shell:
        """
        {params.salmon} quant \
        --gcBias \
        --index {params.index_dir} \
        --libType {params.libtype} \
        --mates1 {input.fast1} \
        --mates2 {input.fast2} \
        --geneMap {input.scallop_ref} \
        --threads {threads} \
        {params.extra_params} \
        -o {params.output_dir} \
        """
