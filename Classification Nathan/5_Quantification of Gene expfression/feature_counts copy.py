rule feature_counts:
	input:
		aligned_bam = star_outdir + "{name}.Aligned.sorted.out.bam",
		aligned_bai = star_outdir + "{name}.Aligned.sorted.out.bam.bai"
	output:
		out_name = feature_counts_outdir + "{name}_featureCounts_results.txt"
	params:
		ref_anno = REFERENCE_ANNOTATION,
		stranded = config['feature_counts_strand_info']
	run:
		if config["end_type"] == "pe":
			shell("{config[feature_counts_path]} -p -t exon -g gene_id -a {params.ref_anno}  --extraAttributes gene_name -o {output.out_name} {params.stranded} {input.aligned_bam}")
		if config["end_type"] == "se":
			shell("{config[feature_counts_path]} -a {params.ref_anno} --extraAttributes gene_name -o {output.out_name} {params.stranded} {input.aligned_bam}")