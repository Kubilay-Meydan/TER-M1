rule RNA_SeqC:
    input:
        gtf = config['reference']['collapsed_gtf'], # dont forget to collaspe the transcript using python script else it wont be proccessed
        bam = config['datadirs']['dedup'] + "/" + "{file}_Aligned.sortedByCoord.out.md.bam"
    output: config['datadirs']['rnaseq_qc'] + "/" + "{file}_Aligned.sortedByCoord.out.md.bam.metrics.tsv",
    params:
        prefix = config['datadirs']['rnaseq_qc'],
    resources:
        mem_mb= 10000
    shell:
        """
         rnaseqc {input.gtf} {input.bam} {params.prefix} --legacy --verbose
        """