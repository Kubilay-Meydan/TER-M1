rule kallisto_quant:
    input:
        config["kallisto_out"] + "/" + config["transcriptome"] + ".idx",
        config["data"] + "/{sample}_1.fastq.gz",
        config["data"] + "/{sample}_2.fastq.gz"
    output:
        config["kallisto_out"] + "/{sample}/abundance.h5",
        config["kallisto_out"] + "/{sample}/abundance.tsv",
        config["kallisto_out"] + "/{sample}/run_info.json"
    conda:
        "../envs/kallisto.yaml"
    params:
        out_dir = config["kallisto_out"] + "/{sample}/bootstraps"
    shell:
        "kallisto quant -i {input[0]}" -o {params.out_dir} -b 100 {input[1]} {input[2]}
    
