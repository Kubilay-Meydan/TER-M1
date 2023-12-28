rule hisat2:
    input:
        hisat2Input
    output:
        sam = pipe('aligned/{sample}.mapped.sam'),
        qc = 'qc/hisat2/{sample}.hisat2.txt'
    params:
        index = config['genome']['index'],
        strand = getHisat2Strand
    log:
        'logs/hisat2/{sample}.log'
    conda:
        f'{ENVS}/hisat2.yaml'
    threads:
        max(1, config['threads'] - 1)
    shell:
        hisat2Cmd()