rule sort_bams:
	input:
		star_outdir + "{name}.Aligned.out.bam"

	output:
		star_outdir + "{name}.Aligned.sorted.out.bam"

	conda:
		"../env/align.yaml"

	shell:
		"""
		samtools sort {input} -o {output}
		"""