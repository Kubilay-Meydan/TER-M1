rule picard_rna_seq:
	input:
		bamfile = expand(star_outdir + "{name}.Aligned.sorted.out.bam",name = SAMPLE_NAMES)
	output:
		config['project_top_level'] + "picard"
	params:
		multiqc_path = config['multiqc_path'],
        rnaseqstrand =
	shell:
		"""
		java -jar /share/apps/genomics/picard-2.20.3/bin/picard.jar CollectRnaSeqMetrics \
		I={input.bamfile} \
		O=output.RNA_Metrics \
		REF_FLAT=ref_flat.txt \
		STRAND= \
		"""