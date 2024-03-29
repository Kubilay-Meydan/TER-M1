rule samtools_flagstat:
	input:
		"bowtie_align/{sample}.sorted.bam"

	output:
		"qc/samtools/{sample}.flagstat"

	threads: 2

	conda :"envs/samtools.yaml"

	shell:"""
		samtools flagstat {input} > {output}
	"""
