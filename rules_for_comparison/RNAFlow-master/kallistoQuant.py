rule kallistoQuant:
    input:
        index = rules.buildKallistoIndex.output,
        reads = rules.cutadapt.output.trimmed
    output:
        qc = 'qc/kallisto/{sample}.stdout',
        abundancesH5 = 'kallisto/{sample}/abundance.h5',
        abundancesTSV = 'kallisto/{sample}/abundance.tsv',
        runInfo = 'kallisto/{sample}/run_info.json'
    params:
        seed = 42,
        fragmentSD = 2,
        bootstraps = 100,
        fragmentLength = 200,
        strand = getKallistoStrand
    log:
        'logs/kallistoQuant/{sample}.log'
    conda:
        f'{ENVS}/kallisto.yaml'
    threads:
        config['threads']
    shell:
        kallistoCmd()