    rule fastqScreen:
        input:
            'fastq/trimmed/{sample}-{read}.trim.fastq.gz'
        output:
            txt = 'qc/fastqScreen/{sample}-{read}.fastq_screen.txt',
            png = 'qc/fastqScreen/{sample}-{read}.fastq_screen.png'
        params:
            config = config['fastqScreen'],
            subset = 100000,
        group:
            'processFASTQ'
        log:
            'logs/fastqScreen/{sample}-{read}.log'
        conda:
            f'{ENVS}/fastqScreen.yaml'
        threads:
            config['threads']
        shell:
            '{SCRIPTS}/fastqScreen.py {input} {params.config} '
            '--subset {params.subset} --threads {threads} '
            '--plotOut {output.png} --dataOut {output.txt} &> {log}'