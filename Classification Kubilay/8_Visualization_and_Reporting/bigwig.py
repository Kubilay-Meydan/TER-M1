## Convert BAM files to bigWig
rule bigwig:
	input:
		bam = os.path.join(outputdir, "STAR", "{sample}", "{sample}_Aligned.sortedByCoord.out.bam"),
		chrl = os.path.join(config["STARindex"], "chrNameLength.txt")
	output:
		os.path.join(outputdir, "STARbigwig", "{sample}_Aligned.sortedByCoord.out.bw")
	params:
		STARbigwigdir = lambda wildcards, output: os.path.dirname(output[0])   ## dirname of first output
	log:
		os.path.join(outputdir, "logs", "bigwig_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "bigwig_{sample}.txt")
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'bedtools version:\n' > {log}; bedtools --version >> {log}; "
		"bedtools genomecov -split -ibam {input.bam} -bg | LC_COLLATE=C sort -k1,1 -k2,2n > "
		"{params.STARbigwigdir}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph; "
		"bedGraphToBigWig {params.STARbigwigdir}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph "
		"{input.chrl} {output}; rm -f {params.STARbigwigdir}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph"