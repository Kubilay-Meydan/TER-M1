rule map:
	input:  sample="raw/{sample}.trimmed.fastq.gz", starref=STARREFDIR, gtf=GTFFILE
	output: "mapped/{sample}.sam"
	threads: 24
	shell:
		"""
		{STAR} --genomeDir {input.starref} --outFileNamePrefix {wildcards.sample}_ --readFilesIn {input.sample} --runThreadN 24 --genomeLoad NoSharedMemory --outSAMattributes All --outSAMstrandField intronMotif --sjdbGTFfile {input.gtf}
		mv {wildcards.sample}_Aligned.out.sam {output}
		mv {wildcards.sample}_Log.final.out {wildcards.sample}_Log.out {wildcards.sample}_Log.progress.out {wildcards.sample}_SJ.out.tab starlogs
		"""
