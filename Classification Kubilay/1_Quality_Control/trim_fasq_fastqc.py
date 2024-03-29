rule trim_fastq_fastqc:
    input:
        pair1 = create_fastq_inputs(config)[0]
    output:
        trimmed_pair1 = temp("output/trim_fastq/{sample}_R1_trimmed.fq.gz"),
        trimmed_pair2 = temp("output/trim_fastq/{sample}_R2_trimmed.fq.gz"),
        fastqc_zipfile1 = "output/fastqc/{sample}_R1_fastqc.zip",
        fastqc_zipfile2 = "output/fastqc/{sample}_R2_fastqc.zip"
    log:
        "output/logs/{sample}.trim_adapters.log"
    params:
        pair2 = create_fastq_inputs(config)[1],
        umi_1 = config["UMI_read1_pattern"],
        umi_2 = config["UMI_read2_pattern"]
    run:
        shell("mkdir -p output/temp_dir")
        if config['trim_polyA'] == "TRUE":
            if config["type"] == "paired":
                if config["use_UMI"] == "TRUE":
                    shell("umi_tools extract -I {input.pair1} --bc-pattern={params.umi_1}  --bc-pattern2={params.umi_2} \
                    --read2-in={params.pair2} --stdout=output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                    --read2-out=output/temp_dir/{wildcards.sample}_R2.fq{suffix}")
                else:
                    # mv files to R1 and R2 ending in temporary directory
                    shell("cp {input.pair1} \
                        output/temp_dir/{wildcards.sample}_R1.fq{suffix}")
                    shell("cp {params.pair2} \
                        output/temp_dir/{wildcards.sample}_R2.fq{suffix}")
                shell("trim_galore \
                    --gzip output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                    output/temp_dir/{wildcards.sample}_R2.fq{suffix} --paired \
                    -o ./output/trim_fastq")
                shell("trim_galore \
                    ./output/trim_fastq/{wildcards.sample}_R1_val_1.fq.gz \
                    ./output/trim_fastq/{wildcards.sample}_R2_val_2.fq.gz --paired --polyA \
                     --basename {wildcards.sample}_pat")
                shell("fastqc output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                    output/temp_dir/{wildcards.sample}_R2.fq{suffix} \
                    -o ./output/fastqc")
                shell("mv ./{wildcards.sample}_pat_val_1.fq.gz \
                    output/trim_fastq/{wildcards.sample}_R1_trimmed.fq.gz"),
                shell("mv ./{wildcards.sample}_pat_val_2.fq.gz \
                    output/trim_fastq/{wildcards.sample}_R2_trimmed.fq.gz")
                shell("mv {wildcards.sample}_R2_val_2.fq.gz_trimming_report.txt ./output/trim_fastq/")
                shell("mv {wildcards.sample}_R1_val_1.fq.gz_trimming_report.txt ./output/trim_fastq/")
                shell("rm ./output/trim_fastq/{wildcards.sample}_R1_val_1.fq.gz")
                shell("rm ./output/trim_fastq/{wildcards.sample}_R2_val_2.fq.gz")
            if config["type"] == "single":
                if config["use_UMI"] == "TRUE":
                    shell("umi_tools extract --stdin={input.pair1} --bc-pattern={params.umi_1} \
                        --log={log} --stdout output/temp_dir/{wildcards.sample}_R1.fq{suffix}")
                else:
                    # mv files to R1 and R2 ending in temporary directory
                    shell("cp {input.pair1} \
                        output/temp_dir/{wildcards.sample}_R1.fq{suffix}")
                shell("trim_galore \
                    --gzip output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                    -o ./output/trim_fastq --basename {wildcards.sample}")
                shell("trim_galore --polyA \
                    ./output/trim_fastq/{wildcards.sample}_trimmed.fq.gz ") # new rule --basename {wildcards.sample}_pat
                shell("mv {wildcards.sample}_trimmed_trimmed.fq.gz \
                    output/trim_fastq/{wildcards.sample}_R1_trimmed.fq.gz")
                shell("mv {wildcards.sample}_trimmed.fq.gz_trimming_report.txt ./output/trim_fastq/")
                shell("fastqc output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                    -o ./output/fastqc")
                shell("touch {output.trimmed_pair2}")
                shell("touch {output.fastqc_zipfile2}")
        else:
            if config["type"] == "paired":
                if config["use_UMI"] == "TRUE":
                    shell("umi_tools extract -I {input.pair1} --bc-pattern={params.umi_1}  --bc-pattern2={params.umi_2} \
                    --read2-in={params.pair2} --stdout=output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                    --read2-out=output/temp_dir/{wildcards.sample}_R2.fq{suffix}")
                else:
                    # mv files to R1 and R2 ending in temporary directory
                    shell("cp {input.pair1} \
                        output/temp_dir/{wildcards.sample}_R1.fq{suffix}")
                    shell("cp {params.pair2} \
                        output/temp_dir/{wildcards.sample}_R2.fq{suffix}")
                if config["experiment"] == "cutrun":
                    shell("trim_galore --clip_R1 6 --clip_R2 6 \
                        --gzip output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                        output/temp_dir/{wildcards.sample}_R2.fq{suffix} --paired --trim-n \
                        -o ./output/trim_fastq")
                else:
                    shell("trim_galore \
                        --gzip output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                        output/temp_dir/{wildcards.sample}_R2.fq{suffix} --paired --trim-n \
                        -o ./output/trim_fastq")
                shell("fastqc output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                    output/temp_dir/{wildcards.sample}_R2.fq{suffix} \
                    -o ./output/fastqc")
                shell("mv output/trim_fastq/{wildcards.sample}_R1_val_1.fq.gz \
                    output/trim_fastq/{wildcards.sample}_R1_trimmed.fq.gz"),
                shell("mv output/trim_fastq/{wildcards.sample}_R2_val_2.fq.gz \
                    output/trim_fastq/{wildcards.sample}_R2_trimmed.fq.gz")
            if config["type"] == "single":
                if config["use_UMI"] == "TRUE":
                    shell("umi_tools extract --stdin={input.pair1} --bc-pattern={params.umi_1} \
                        --log={log} --stdout output/temp_dir/{wildcards.sample}_R1.fq{suffix}")
                else:
                    # mv files to R1 and R2 ending in temporary directory
                    shell("cp {input.pair1} \
                        output/temp_dir/{wildcards.sample}_R1.fq{suffix}")
                shell("trim_galore --trim-n \
                    --gzip output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                    -o ./output/trim_fastq --basename {wildcards.sample}")
                shell("mv output/trim_fastq/{wildcards.sample}_trimmed.fq.gz \
                    output/trim_fastq/{wildcards.sample}_R1_trimmed.fq.gz")
                shell("fastqc output/temp_dir/{wildcards.sample}_R1.fq{suffix} \
                    -o ./output/fastqc")
                shell("touch {output.trimmed_pair2}")
                shell("touch {output.fastqc_zipfile2}")
