rule featureCounts:
    input:
        bam = rules.markdupBAM.output.bam,
        gtf = rules.gff3ToGTF.output
    output:
        counts = 'featureCounts/{sample}.featureCounts.txt',
        qc = 'featureCounts/{sample}.featureCounts.txt.summary'
    params:
        paired = '-pC' if config['paired'] else '',
        strand = getStrand
    log:
        'logs/featureCounts/{sample}.log'
    conda:
        f'{ENVS}/subread.yaml'
    threads:
        config['threads']
    shell:
        'featureCounts {params.paired} -a {input.gtf} -o {output.counts} '
        '{input.bam} -t exon -g gene_id -T {threads} -s {params.strand} '
        '&> {log}'

