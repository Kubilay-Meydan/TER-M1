rule fastQCtrimmed:
    input:
        'fastq/trimmed/{sample}-{read}.trim.fastq.gz'
    output:
        html = 'qc/fastqc/{sample}-{read}.trim_fastqc.html',
        zip = 'qc/fastqc/{sample}-{read}.trim_fastqc.zip'
    group:
        'processFASTQ'
    log:
        'logs/fastqc_trimmed/{sample}-{read}.log'
    wrapper:
        '0.49.0/bio/fastqc'