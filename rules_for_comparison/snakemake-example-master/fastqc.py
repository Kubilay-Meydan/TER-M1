rule fastqc:
	input: "raw/{sample}.trimmed.fastq.gz"
	output: "fastqc/{sample}.trimmed_fastqc.zip"
	shell: "{TOOLDIR}/FastQC/fastqc -o fastqc {input}"
	