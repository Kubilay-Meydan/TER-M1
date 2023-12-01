rule markdupBAM:
    input:
        'aligned/{sample}.sort.bam'
    output:
        bam = 'aligned/{sample}.markdup.bam',
        qc = 'qc/markDuplicates/{sample}.markDuplicates.log'
    params:
        tmp = config['tmpdir']
    log:
        'logs/markDuplicates/{sample}.log'
    conda:
        f'{ENVS}/picard.yaml'
    shell:
        'picard MarkDuplicates I={input} O={output.bam} M={output.qc} '
        'VALIDATION_STRINGENCY=LENIENT TMP_DIR={params.tmp} &> {log}'
