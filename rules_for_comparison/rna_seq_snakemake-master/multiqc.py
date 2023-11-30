rule multiqc:
    input:
        multiqc_target_files(workflow_str, SAMPLE_NAMES, FASTQ_PREFIX, UNITS)

    output:
        os.path.join(multiqc_output_folder, "multiqc_report.html")

    params:
        dirs = multiqc_target_dirs(),
        outdir = multiqc_output_folder,
        save_plots = "-p", # Should provide option to alter this in config
        runtime = "--profile-runtime"

    conda:
        "../env/align.yaml"

    shell:
        """
        multiqc \
        {params.save_plots} \
        -o {params.outdir} \
        {params.runtime} \
        {params.dirs}
        """