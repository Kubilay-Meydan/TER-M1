rule featureCounts:
	input:
		bam="filtered_bam/{sample}.bam",
		gff="refs/mirbase.gff",
		fasta="refs/reference.fasta"

	output:
		multiext("featureCounts/{sample}",".featureCounts", ".featureCounts.summary")

	threads: 2

	conda : "envs/subread.yaml"

	log: "featureCounts/{sample}.txt"

	shell:"""

		featureCounts \
		-T {threads} \
		-t miRNA \
		-g Name \
		-G {input.fasta} \
		-F GFF \
		--fracOverlap 0.2 \
		-O \
		-s 0 \
		-M \
		-a \
		{input.gff} \
		-o {output[0]} {input.bam} 2> {log}
	"""
