rule salmon_quant:
	input:
		fastq="filtered_fastq/{sample}.fastq.gz",
		index=directory("refs/salmon_index")

	output:
		quant="salmon/{sample}_quant/quant.sf",

	threads: 4

	params:
		outdir="Salmon/{sample}"

	conda: "envs/salmon.yaml"

	shell:"""
		salmon quant \
		-i {input.index} \
		-l A \
		-p {threads} \
		-r {input.fastq} \
		--validateMappings \
		-o {params.outdir}_quant
	"""