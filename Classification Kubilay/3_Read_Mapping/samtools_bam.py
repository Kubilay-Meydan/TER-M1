rule samtools_bam:
	input:
		"bowtie_align/{sample}.sam"

	output:
		temp("bowtie_align/{sample}.bam")

	threads: 2

	conda : "envs/samtools.yaml"

	shell:"""
		samtools view -bS {input} > {output}

	"""
