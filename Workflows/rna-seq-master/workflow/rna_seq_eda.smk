rule make_salmon_txi:
    input: expand(salmon_dir + "/{library}.quant.sf", library = RNA_LIBS),
    log: logdir + "/{experiment}_make_salmon_txi.log",
    output: analysis + "/{experiment}_txi.rdata",
    params:
        script = rna_script_dir + "/make_salmon_txi.R",
        txdb = txdb,
    shell:
        """
        Rscript {params.script} \
        "{input}" \
        {output} \
        {params.txdb} \
        > {log} 2>&1
        """

rule all_rna_eda:
    container: "/home/jeszyman/sing_containers/atac.1.1.0.sif",
    input: analysis + "/{experiment}_txi.rdata",
    log: logdir + "/{experiment}_rna_eda.log",
    output:
        pca = results + "/figures/{experiment}_pca.pdf",
        rdata = analysis + "/{experiment}_eda.rdata",
    params:
        factor_str = factor_str,
        library_tsv = library_tsv,
        script = rna_scriptdir + "/all_rna_eda.R",
    shell:
        """
        Rscript {params.script} \
        {input} \
        {output.pca} \
        {output.rdata} \
        "{params.factor_str}" \
        {params.library_tsv} \
        > {log} 2>&1
        """
