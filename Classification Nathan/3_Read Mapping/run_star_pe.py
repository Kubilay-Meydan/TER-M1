rule run_star_pe:
	wildcard_constraints:
		sample="|".join(SAMPLE_NAMES)
	input:
		generated_index = GENOME_DIR + "/SA",
		one = lambda wildcards: get_processed_fastq(wildcards.name, pair=1),
		two = lambda wildcards: get_processed_fastq(wildcards.name, pair=2)
	output:
		star_outdir + "{name}.SJ.out.tab",
		star_outdir + "{name}.Log.final.out",
		temp(star_outdir + "{name}.Aligned.out.bam")
	params:
		extra_star_parameters = return_parsed_extra_params(config['extra_star_parameters']),
		genomeDir = GENOME_DIR,
		outTmpDir = os.path.join(star_outdir + "{name}_tmpdir"),
		outputPrefix = os.path.join(star_outdir + "{name}.")
	threads:
		4
	conda:
		"../env/align.yaml"

	shell:
		"""
		rm -rf {params.outTmpDir}
		STAR --genomeDir {params.genomeDir} \
		--readFilesIn {input.one} {input.two} \
		--outFileNamePrefix {params.outputPrefix} \
		--readFilesCommand zcat --runThreadN {threads} \
		{params.extra_star_parameters} \
		--outTmpDir {params.outTmpDir}
		"""