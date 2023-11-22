rule trimgalorePE:
	input:
		fastq1 = os.path.join(FASTQdir, "".join(["{sample}_", str(config["fqext1"]), ".", str(config["fqsuffix"]), ".gz"])),
		fastq2 = os.path.join(FASTQdir, "".join(["{sample}_", str(config["fqext2"]), ".", str(config["fqsuffix"]), ".gz"]))
	output:
		os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext1"]), "_val_1.fq.gz"])),
		os.path.join(outputdir, "FASTQtrimmed", "".join(["{sample}_", str(config["fqext2"]), "_val_2.fq.gz"]))
	params:
		FASTQtrimmeddir = lambda wildcards, output: os.path.dirname(output[0])   ## dirname of first output
	log:
		os.path.join(outputdir, "logs", "trimgalore_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "trimgalore_{sample}.txt")
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'TrimGalore! version:\n' > {log}; trim_galore --version >> {log}; "
		"trim_galore -q 20 --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt "
		"--paired {input.fastq1} {input.fastq2}"