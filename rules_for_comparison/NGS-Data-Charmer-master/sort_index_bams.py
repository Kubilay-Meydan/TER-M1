rule sort_index_bams:
	input:
		star_outdir + "{name}.Aligned.sorted.out.bam"

	output:
		star_outdir + "{name}.Aligned.sorted.out.bam.bai"

	conda:
		"../env/align.yaml"

	shell:
		"""
		samtools index {input}
		"""