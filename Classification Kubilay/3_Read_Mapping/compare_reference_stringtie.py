rule compare_reference_stringtie:
    input:
        os.path.join(stringtie_outdir,"stringtie_merged.gtf")
    output:
        os.path.join(stringtie_outdir, "gffall.stringtie_merged.gtf.tmap")
    params:
        ref_gtf = GTF,
        gffcompare = config['gffcompare']
    shell:
        """
        {params.gffcompare} -o gffall -r {params.ref_gtf} {input}
        """
