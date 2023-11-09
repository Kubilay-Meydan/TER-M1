rule kallisto_index:
    input:
        config["transcriptome"]
    output:
        config["kallisto_out"] + "/" + config["transcriptome"] + ".idx"
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto -i {output[0]} {input[0]}"
