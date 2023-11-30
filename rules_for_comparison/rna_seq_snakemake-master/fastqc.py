rule fastqc:
	input:
		fastq_file = lambda wildcards: return_fastq_location(wildcards.sample),
		#fq_name = lambda wildcards, input: re.sub(".fastq.gz","", ORDER_DICT[input.fastq_file].rpartition('/')[2]),
	output:
		out_fastqc = fastqc_outdir + "{sample}/{unit}/{fastq_prefix}_fastqc.html"
	params:
		#fq_name = lambda wildcards, input: re.sub(".fastq.gz","", ORDER_DICT[input.fastq_file].rpartition('/')[2]),
		outdir = fastqc_outdir + "{sample}/{unit}/"

	conda:
		"../env/align.yaml"

	shell:
		"""
		mkdir -p {params.outdir}
		fastqc {input.fastq_file} -o {params.outdir}
		"""
