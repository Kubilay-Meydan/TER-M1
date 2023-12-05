
__author__ = "Masood Zaka (https://github.com/masoodzaka/Snakemake_smallRNASeq.git)"
__licence__ = "MIT"

shell.prefix("set -o pipefail; ")

from os import path
import re
from snakemake.utils import validate
from snakemake.utils import min_version

SAMPLES=["DFU01","DFU02","DFU03","DFU04","DM01","DM02","DM03","DM04"]

rule all:
	input:"qc/multiqc_report.html",
		expand("salmon/{sample}_quant/quant.sf", sample=SAMPLES)


rule get_genome:
	output:
		dna="refs/reference.fasta",
		trans="refs/transcripts.fasta",
		gentrome="refs/gentrome.fasta"

	shell:"""
			curl -L ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz | gzip -d > {output.dna}
			curl -L ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.transcripts.fa.gz | gzip -d > {output.trans}
			cat {output.trans} {output.dna} > {output.gentrome}

		"""
rule get_annotation:
	output:
		"refs/reference.gtf"

	shell:"""
		curl -L ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz | gzip -d > {output}
	"""
rule get_gff:
	output:
		gff="refs/mirbase.gff",
		bed="refs/mirbase.bed"
	shell:"""
		curl -L https://www.mirbase.org/ftp/22/genomes/hsa.gff3 | awk '{{if ($3=="miRNA") print}}' > {output.gff}
		curl -L https://www.mirbase.org/ftp/22/genomes/hsa.gff3 | awk 'BEGIN {{OFS = "\t"}};{{if ($3=="miRNA") print$1,$4-1,$5}}' > {output.bed}
	"""

rule fastqc:
    input:
        "data/Reads/{sample}/{sample}.fq.gz"
    output:
        html="qc/fastqc/{sample}_fastqc.html",
        zip="qc/fastqc/{sample}_fastqc.zip"

    threads: 2

    params:
        "--quiet"

    conda: "envs/qc.yaml",

    log: "qc/fastqc/{sample}.log"

	shell:"""
        fastqc -t {threads} {input} --noextract {params} -o qc/fastqc 2> {log}
        """

rule cutadapt:
	input:
		reads=("data/Reads/{sample}/{sample}.fq.gz"),

	output:
		"cutadapt/{sample}.trimmed.fq.gz"

	threads: 2

	log: "qc/cutadapt/{sample}.cutadapt.json"

	conda: "envs/cutadapt.yaml"

	shell:"""
		cutadapt \
		--cores {threads} \
		-q 20 \
		-m 18 \
		-M 30 \
		--json {log} \
		-o {output} \
		{input}
	"""
rule fastqc_trimmed:
    input:
        "cutadapt/{sample}.trimmed.fq.gz"
    output:
        html="qc/fastqc_trimmed/{sample}.trimmed_fastqc.html",
        zip="qc/fastqc_trimmed/{sample}.trimmed_fastqc.zip"

    threads: 2

    params:
        "--quiet"

    conda: "envs/qc.yaml",

    log: "qc/fastqc_trimmed/{sample}.log"

	shell:"""
        fastqc -t {threads} {input} --noextract {params} -o qc/fastqc_trimmed 2> {log}
        """

rule bowtie_build:
	input:
		"refs/reference.fasta"

	output:
		index=multiext("refs/reference",".1.ebwt",".2.ebwt",".3.ebwt",".4.ebwt",".rev.1.ebwt",".rev.2.ebwt")

	threads: 8

	params:
		index=lambda w, input: os.path.splitext(input[0])[0]

	conda: "envs/bowtie.yaml"

	shell:"""
		bowtie-build \
		--threads {threads} \
		{input} \
		{params.index}
	"""
rule bowtie_align:
	input:
		reads=("cutadapt/{sample}.trimmed.fq.gz"),
		index=multiext("refs/reference",".1.ebwt",".2.ebwt",".3.ebwt",".4.ebwt",".rev.1.ebwt",".rev.2.ebwt")

	output:
		temp("bowtie_align/{sample}.sam")

	threads: 8

	params:
		index=lambda w, input: os.path.commonprefix(input.index).rstrip(".")

	conda:"envs/bowtie.yaml"

	log: "bowtie_align/{sample}.log"

	shell:"""
		bowtie \
		--threads {threads} \
		-n 1 \
		-l 20 \
		--best \
		-m 1 \
		--sam \
		{input.reads} \
		-x {params.index} \
		{output} 2>{log}
	"""

rule samtools_bam:
	input:
		"bowtie_align/{sample}.sam"

	output:
		temp("bowtie_align/{sample}.bam")

	threads: 2

	conda : "envs/samtools.yaml"

	shell:"""
		samtools view -bS {input} > {output}

	"""

rule samtools_sort:
	input:
		"bowtie_align/{sample}.bam"

	output:
		"bowtie_align/{sample}.sorted.bam"

	threads: 2

	conda : "envs/samtools.yaml"

	shell:"""
		samtools sort -o {output} {input}

	"""
rule samtools_index:
	input:
		"bowtie_align/{sample}.sorted.bam"

	output:
		"bowtie_align/{sample}.sorted.bam.bai"

	threads: 2

	conda : "envs/samtools.yaml"

	shell:"""
		samtools index {input}

	"""
rule samtools_flagstat:
	input:
		"bowtie_align/{sample}.sorted.bam"

	output:
		"qc/samtools/{sample}.flagstat"

	threads: 2

	conda :"envs/samtools.yaml"

	shell:"""
		samtools flagstat {input} > {output}
	"""

rule samtools_filter_bam:
	input:
		bam="bowtie_align/{sample}.sorted.bam",
		bed="refs/mirbase.bed"

	output:
		"filtered_bam/{sample}.bam"

	threads: 2

	conda :"envs/samtools.yaml"

	shell:"""
		samtools view -b -h -L {input.bed} {input.bam}> {output}
	"""
rule samtools_filter_bam_index:
	input:
		"filtered_bam/{sample}.bam"

	output:
		"filtered_bam/{sample}.bam.bai"

	threads: 2

	conda :"envs/samtools.yaml"

	shell:"""
		samtools index {input}
	"""
rule bedtools_annotate:
	input:
		bam="filtered_bam/{sample}.bam",
		gff="refs/mirbase.gff",

	output:
		"annotation/{sample}.bed"

	threads: 1

	conda : "envs/bedtools.yaml"

	shell:"""

		intersectBed -abam {input.bam} -b {input.gff} -r -wa -wb -bed > {output}
	"""

rule featureCounts:
	input:
		bam="filtered_bam/{sample}.bam",
		gff="refs/mirbase.gff",
		fasta="refs/reference.fasta"

	output:
		multiext("featureCounts/{sample}",".featureCounts", ".featureCounts.summary")

	threads: 2

	conda : "envs/subread.yaml"

	log: "featureCounts/{sample}.txt"

	shell:"""

		featureCounts \
		-T {threads} \
		-t miRNA \
		-g Name \
		-G {input.fasta} \
		-F GFF \
		--fracOverlap 0.2 \
		-O \
		-s 0 \
		-M \
		-a \
		{input.gff} \
		-o {output[0]} {input.bam} 2> {log}
	"""

rule samtools_bam2fastq:
	input:
		"filtered_bam/{sample}.bam"

	output:
		"filtered_fastq/{sample}.fastq.gz"

	threads: 2

	conda :"envs/samtools.yaml"

	shell:"""
		samtools fastq {input} | gzip > {output}
	"""

rule salmon_Index:
    input:
        "refs/transcripts.fasta",

    output:
        index=directory("refs/salmon_index"),

    threads: 8

    conda: "envs/salmon.yaml"

    shell:""" salmon \
            index \
            -k 21 \
            --threads {threads} \
            -t {input} \
            --gencode \
            -i {output.index}

"""
rule salmon_quant:
	input:
		fastq="filtered_fastq/{sample}.fastq.gz",
		index=directory("refs/salmon_index")

	output:
		quant="salmon/{sample}_quant/quant.sf",

	threads: 4

	params:
		outdir="Salmon/{sample}"

	conda: "envs/salmon.yaml"

	shell:"""
		salmon quant \
		-i {input.index} \
		-l A \
		-p {threads} \
		-r {input.fastq} \
		--validateMappings \
		-o {params.outdir}_quant
	"""
rule multiqc:
    input:
	    expand(["qc/fastqc/{sample}_fastqc.zip",
	    "qc/cutadapt/{sample}.cutadapt.json",
	    "qc/fastqc_trimmed/{sample}.trimmed_fastqc.zip",
	    "bowtie_align/{sample}.log",
		"qc/samtools/{sample}.flagstat",
		"featureCounts/{sample}.featureCounts.summary",
		"salmon/{sample}_quant/"], sample=SAMPLES)

    output:
        "qc/multiqc_report.html"

    threads: 2

    conda: "envs/qc.yaml"

    shell:"""
        multiqc {input} -n multiqc_report -f -q -o qc
    """

