rule multiqc:
	input:
		multiqc_input
	output:
		os.path.join(outputdir, "MultiQC", "multiqc_report.html")
	params:
		inputdirs = multiqc_params,
		MultiQCdir = lambda wildcards, output: os.path.dirname(output[0])   ## dirname of first output
	log:
		os.path.join(outputdir, "logs", "multiqc.log")
	benchmark:
		os.path.join(outputdir, "benchmarks", "multiqc.txt")
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'MultiQC version:\n' > {log}; multiqc --version >> {log}; "
		"multiqc {params.inputdirs} -f -o {params.MultiQCdir}"