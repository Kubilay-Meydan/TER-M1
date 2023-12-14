rule salmonPE:
	input:
		index = os.path.join(config["salmonindex"], "versionInfo.json"),
		fastq1 = os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext1"]), "_val_1.fq.gz"])) if config["run_trimming"] else os.path.join(FASTQdir, "".join(["{sample}_", str(config["fqext1"]), ".", str(config["fqsuffix"]), ".gz"])),
		fastq2 = os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext2"]), "_val_2.fq.gz"])) if config["run_trimming"] else os.path.join(FASTQdir, "".join(["{sample}_", str(config["fqext2"]), ".", str(config["fqsuffix"]), ".gz"]))
	output:
		os.path.join(outputdir, "salmon", "{sample}", "quant.sf")
	log:
		os.path.join(outputdir, "logs", "salmon_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "salmon_{sample}.txt")
	threads:
		config["ncores"]
	params:
		salmonindex = lambda wildcards, input: os.path.dirname(input['index']),   ## dirname of index input
		salmondir = lambda wildcards, output: os.path.dirname(os.path.dirname(output[0])),   ## dirname of first output
		salmonextraparams = config["additional_salmon_quant"]
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'Salmon version:\n' > {log}; salmon --version >> {log}; "
		"salmon quant -i {params.salmonindex} -l A -1 {input.fastq1} -2 {input.fastq2} "
		"-o {params.salmondir}/{wildcards.sample} -p {threads} {params.salmonextraparams}"
