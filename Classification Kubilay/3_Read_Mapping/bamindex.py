## Index bam files
rule bamindex:
	input:
		bam = os.path.join(outputdir, "STAR", "{sample}", "{sample}_Aligned.sortedByCoord.out.bam")
	output:
		os.path.join(outputdir, "STAR", "{sample}", "{sample}_Aligned.sortedByCoord.out.bam.bai")
	log:
		os.path.join(outputdir, "logs", "samtools_index_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "samtools_index_{sample}.txt")
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'samtools version:\n' > {log}; samtools --version >> {log}; "
		"samtools index {input.bam}"