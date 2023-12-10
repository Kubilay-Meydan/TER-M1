rule starPE:
	input:
		index = os.path.join(config["STARindex"], "SA"),
		fastq1 = os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext1"]), "_val_1.fq.gz"])) if config["run_trimming"] else os.path.join(FASTQdir, "".join(["{sample}_", str(config["fqext1"]), ".", str(config["fqsuffix"]), ".gz"])),
		fastq2 = os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext2"]), "_val_2.fq.gz"])) if config["run_trimming"] else os.path.join(FASTQdir, "".join(["{sample}_", str(config["fqext2"]), ".", str(config["fqsuffix"]), ".gz"]))
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
		"STAR --genomeDir {params.STARindex} --readFilesIn {input.fastq1} {input.fastq2} "
		"--runThreadN {threads} --outFileNamePrefix {params.STARdir}/{wildcards.sample}/{wildcards.sample}_ "
		"--outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c "
		"{params.starextraparams}"