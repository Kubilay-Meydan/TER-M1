rule fastqc_trimmed:
    input:
        "cutadapt/{sample}.trimmed.fq.gz"
    output:
        html="qc/fastqc_trimmed/{sample}.trimmed_fastqc.html",
        zip="qc/fastqc_trimmed/{sample}.trimmed_fastqc.zip"

    threads: 2

    params:
        "--quiet"

    conda: "envs/qc.yaml",

    log: "qc/fastqc_trimmed/{sample}.log"

	shell:"""
        fastqc -t {threads} {input} --noextract {params} -o qc/fastqc_trimmed 2> {log}
        """
