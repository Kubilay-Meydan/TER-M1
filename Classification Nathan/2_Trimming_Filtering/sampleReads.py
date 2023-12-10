rule sampleReads:
    input:
        'fastq/trimmed/{sample}-{read}.trim.fastq.gz'
    output:
        'fastq/sampled/{sample}-{read}.trim.fastq.gz'
    params:
        seed = 42,
        nReads = config['misc']['sample']
    log:
        'logs/sampleReads/{sample}-{read}.log'
    conda:
        f'{ENVS}/seqtk.yaml'
    shell:
        '(seqtk sample -s {params.seed} {input} {params.nReads} '
        '| gzip > {output}) 2> {log}'