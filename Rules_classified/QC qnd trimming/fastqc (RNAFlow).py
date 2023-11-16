rule fastqc:
    input:
        lambda wc: RNASeq.path(wc.sample, [wc.read])
    output:
        html = 'qc/fastqc/{sample}-{read}.raw_fastqc.html',
        zip = 'qc/fastqc/unmod/{sample}-{read}.raw.fastqc.zip'
    group:
        'processFASTQ'
    log:
        'logs/fastqc/{sample}-{read}.log'
    wrapper:
        '0.49.0/bio/fastqc'
