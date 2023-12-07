rule salmon_index_extended:
    input:
        extended_fa = os.path.join(scallop_outdir, "gentrome.fa"),
        decoys = os.path.join(scallop_outdir, "decoys.txt")
    output:
        os.path.join(scallop_outdir, "extended_transcriptome/seq.bin"),
        os.path.join(scallop_outdir, "extended_transcriptome/pos.bin")
    params:
        salmon = config["salmon_path"],
        k = KMER_SIZE,
        outdir = os.path.join(scallop_outdir, "extended_transcriptome"),
        gencode = "--gencode" if config["transcriptome_source"] == "gencode" else ""
    threads:
        4
    shell:
        """
        {params.salmon} index \
        -t {input.extended_fa} \
        -i {params.outdir} \
        --decoys {input.decoys} \
        -k {params.k} \
        {params.gencode} \
        -p {threads}
        """