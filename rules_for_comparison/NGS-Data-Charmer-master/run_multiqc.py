
# Create multiqc report (used for all workflows) 
rule run_multiqc:
    input:
        "output/fpkm_genic_matrix.txt" if config["experiment"] == "rnaseq" else \
        expand("output/bam/{sample}.unique.sorted.rmdup.chr.bam", sample=SAMPLES)
    output:
        multiqc_report = "output/multiqc_report.html"
    params:
        multiqc_config = os.path.join(expand("{param}", param=config["ngs_path"])[0],"multiqc_config_template.yaml")
    shell:
        "multiqc . -f --outdir ./ --config {params.multiqc_config}"

