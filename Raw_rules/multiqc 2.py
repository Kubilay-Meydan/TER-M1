rule multiqc:
    input:
	    expand(["qc/fastqc/{sample}_fastqc.zip",
	    "qc/cutadapt/{sample}.cutadapt.json",
	    "qc/fastqc_trimmed/{sample}.trimmed_fastqc.zip",
	    "bowtie_align/{sample}.log",
		"qc/samtools/{sample}.flagstat",
		"featureCounts/{sample}.featureCounts.summary",
		"salmon/{sample}_quant/"], sample=SAMPLES)

    output:
        "qc/multiqc_report.html"

    threads: 2

    conda: "envs/qc.yaml"

    shell:"""
        multiqc {input} -n multiqc_report -f -q -o qc
    """
