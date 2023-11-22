rule trimgaloreSE:
	input:
		fastq = os.path.join(FASTQdir, "".join(["{sample}.", str(config["fqsuffix"]), ".gz"]))
	output:
		os.path.join(outputdir, "FASTQtrimmed", "{sample}_trimmed.fq.gz")
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
		"trim_galore -q 20 --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt {input.fastq}"