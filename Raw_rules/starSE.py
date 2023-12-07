## Genome mapping with STAR
rule starSE:
	input:
		index = os.path.join(config["STARindex"], "SA"),
		fastq = os.path.join(outputdir, "FASTQtrimmed", "{sample}_trimmed.fq.gz") if config["run_trimming"] else os.path.join(FASTQdir, "".join(["{sample}.", str(config["fqsuffix"]), ".gz"]))
	output:
		os.path.join(outputdir, "STAR", "{sample}", "{sample}_Aligned.sortedByCoord.out.bam")
	threads:
		config["ncores"]
	log:
		os.path.join(outputdir, "logs", "STAR_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "STAR_{sample}.txt")
	params:
		STARindex = lambda wildcards, input: os.path.dirname(input['index']),   ## dirname of index input
		STARdir = lambda wildcards, output: os.path.dirname(os.path.dirname(output[0])),   ## dirname of first output
		starextraparams = config["additional_star_align"]
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'STAR version:\n' > {log}; STAR --version >> {log}; "
		"STAR --genomeDir {params.STARindex} --readFilesIn {input.fastq} "
		"--runThreadN {threads} --outFileNamePrefix {params.STARdir}/{wildcards.sample}/{wildcards.sample}_ "
		"--outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c "
		"{params.starextraparams}"
