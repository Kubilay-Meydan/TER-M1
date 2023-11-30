rule fetch_unique_stringtie:
    input:
        sample_tmap = os.path.join(stringtie_outdir,"stringtie_merged.gtf"),
        sample_gtf = os.path.join(stringtie_outdir, "gffall.stringtie_merged.gtf.tmap")
    output:
        os.path.join(stringtie_outdir, "stringtie_merged.unique.gtf")
    params:
        ref_gtf = GTF,
        gtfcuff = config['gtfcuff']
    shell:
        """
        {params.gtfcuff} puniq {input.sample_tmap} {input.sample_gtf} {params.ref_gtf} {output}
        """