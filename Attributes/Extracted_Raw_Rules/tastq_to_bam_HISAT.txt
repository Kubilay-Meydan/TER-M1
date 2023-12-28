rule fastq_to_bam_HISAT:
    input:
        trimmed_pair = ["output/trim_fastq/{sample}_R1_trimmed.fq.gz", \
                            "output/trim_fastq/{sample}_R2_trimmed.fq.gz"]
    params:
        hisat2index = config["hisat2_index"]
    output:
        bam = "output/bam/{sample}.bam",
        bambai = "output/bam/{sample}.bam.bai"
    threads: config["threads_for_alignment"]
    log:
        "output/logs/{sample}.alignment.log"
    run:
        timelog="output/logs/"+wildcards.sample+".alignment.log"
        # print(timelog)
        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
        starttime = "Alignment starting date and time =" + dt_string
        ## Timestamp in log file, using datetime object containing current date and time
        ## ChIP-Seq alignment rules
        if config["experiment"] != "rnaseq" :
            if config["type"] == "paired":
                # Perform the alignment
                # Splicing is not desired in cut&run and chipseq
                if config["cufflinks_bam"] == "FALSE" :
                    shell("hisat2 -p {threads} {sensitivity_level} -x {params.hisat2index} \
                            -1 {input.trimmed_pair[0]} \
                            -2 {input.trimmed_pair[1]} \
                            --no-spliced-alignment \
                            -S output/bam/{wildcards.sample}.sam 2> {log}")
                else:
                    shell("hisat2 -p {threads} {sensitivity_level} -x {params.hisat2index} \
                            --pen-noncansplice 1000000 \
                            --no-spliced-alignment \
                            -1 {input.trimmed_pair[0]} \
                            -2 {input.trimmed_pair[1]} \
                            -S output/bam/{wildcards.sample}.sam 2> {log}")
            if config["type"] == "single":        
                if config["cufflinks_bam"] == "FALSE":
                    shell("hisat2 -p {threads} {sensitivity_level} -x {params.hisat2index} \
                            -U {input.trimmed_pair[0]} \
                            --no-spliced-alignment \
                            -S output/bam/{wildcards.sample}.sam 2> {log}")
                else:
                    shell("hisat2 -p {threads} {sensitivity_level} -x {params.hisat2index} \
                            --pen-noncansplice 1000000 -U {input.trimmed_pair[0]} \
                            --no-spliced-alignment \
                            -S output/bam/{wildcards.sample}.sam 2> {log}") 
        ## RNA-seq alignment rules
        if config["experiment"] == "rnaseq" :
            # Splicing is desired in rnaseq
            if config["type"] == "paired":
                if config["cufflinks_bam"] == "FALSE" :
                    shell("hisat2 -p {threads} {sensitivity_level} -x {params.hisat2index} \
                            -1 {input.trimmed_pair[0]} -2 {input.trimmed_pair[1]} \
                            -S output/bam/{wildcards.sample}.sam 2> {log}")
                else:
                    shell("hisat2 -p {threads} {sensitivity_level} -x {params.hisat2index} \
                            --pen-noncansplice 1000000 \
                            -1 {input.trimmed_pair[0]} -2 {input.trimmed_pair[1]} \
                            -S output/bam/{wildcards.sample}.sam 2> {log}")
            if config["type"] == "single":        
                if config["cufflinks_bam"] == "FALSE":
                    shell("hisat2 -p {threads} {sensitivity_level} -x {params.hisat2index} \
                            -U {input.trimmed_pair[0]} \
                            -S output/bam/{wildcards.sample}.sam 2> {log}")
                else:
                    shell("hisat2 -p {threads} {sensitivity_level} -x {params.hisat2index} \
                            --pen-noncansplice 1000000 -U {input.trimmed_pair[0]} \
                            -S output/bam/{wildcards.sample}.sam 2> {log}")
        # Timestamp in log file
        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
        endtime = "Alignment ending date and time =" + dt_string
        with open(timelog, 'a') as f:
            print(starttime, file=f)
            print(endtime, file=f)
        #
        ## Convert and cleanup the alignment files                             
        shell("samtools sort -@ 8 -O BAM -o {output.bam} output/bam/{wildcards.sample}.sam")
        shell("rm output/bam/{wildcards.sample}.sam")
        shell("samtools index {output.bam}")
        if config["experiment"] != "cutrun" :
            shell("rm output/temp_dir/{wildcards.sample}_R1.fq{suffix}")
            if config["type"] == "paired":
                shell("rm output/temp_dir/{wildcards.sample}_R2.fq{suffix}")
        if config["keep_fastq"] == "FALSE":
            shell("rm output/trim_fastq/{wildcards.sample}_R1_trimmed.fq.gz \
                    output/trim_fastq/{wildcards.sample}_R2_trimmed.fq.gz")
