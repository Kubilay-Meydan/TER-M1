rule expressbam:
	input: fq="raw/{sample}.trimmed.fastq.gz", ref=CDNA+'.nix'
	output: "express/bams/{sample}.sam"
	log: "logs/express/{sample}.log"
	shell: "{ALIGN} -d {input.ref} -rALL -f {input.fq} -o SAM 2> {log} > {output}"
