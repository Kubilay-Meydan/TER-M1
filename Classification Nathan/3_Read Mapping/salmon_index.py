rule salmon_index:
    input:
        gentrome_fa = os.path.join(DECOYS_DIR, "gentrome.fa"),
        decoys = os.path.join(DECOYS_DIR, "decoys.txt")

    output:
        os.path.join(TXOME_DIR, "seq.bin"),
        os.path.join(TXOME_DIR, "pos.bin")

    params:
        k = KMER_SIZE,
        outdir = TXOME_DIR,
        gencode = "--gencode" if config["transcriptome_source"] == "gencode" else ""

    threads:
        4

    conda:
        "../env/align.yaml"

    shell:
        """
        salmon index \
        -t {input.gentrome_fa} \
        -i {params.outdir} \
        --decoys {input.decoys} \
        -k {params.k} \
        {params.gencode} \
        -p {threads}
        """
