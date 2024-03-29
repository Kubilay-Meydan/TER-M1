## FastQC, trimmed reads
rule fastqctrimmed:
	input:
		fastq = os.path.join(outputdir, "FASTQtrimmed", "{sample}.fq.gz")
	output:
		os.path.join(outputdir, "FastQC", "{sample}_fastqc.zip")
	params:
		FastQC = lambda wildcards, output: os.path.dirname(output[0])   ## dirname of first output
	log:
		os.path.join(outputdir, "logs", "fastqc_trimmed_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "fastqc_trimmed_{sample}.txt")
	conda:
		"envs/environment.yaml"
	threads:
		config["ncores"]
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		"fastqc -o {params.FastQC} -t {threads} {input.fastq}"
