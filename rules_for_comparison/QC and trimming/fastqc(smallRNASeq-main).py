rule fastqc:
    input:
        "data/Reads/{sample}/{sample}.fq.gz"
    output:
        html="qc/fastqc/{sample}_fastqc.html",
        zip="qc/fastqc/{sample}_fastqc.zip"

    threads: 2

    params:
        "--quiet"

    conda: "envs/qc.yaml",

    log: "qc/fastqc/{sample}.log"

	shell:"""
        fastqc -t {threads} {input} --noextract {params} -o qc/fastqc 2> {log}
        """