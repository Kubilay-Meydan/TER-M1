rule sortedbam_to_rmdup:
        input:
            "output/bam/{sample}.bam"
        output:
            "output/bam/{sample}.sorted.bam" if config["type"] == "single" else "output/bam/{sample}.sorted.rmdup.bam"
        log:
            "output/logs/{sample}.rmdup.log"
        params:
            refflat = config["refflat"]
        run:
            if config["type"] == "paired":
                if config["use_UMI"] == "TRUE":
                    shell("umi_tools dedup --stdin={input} --log={log} --paired --unmapped-reads=use > {output}")
                else:
                    # shell("samtools markdup -r {input} {output} 2> {log}")
                    shell("picard MarkDuplicates REMOVE_DUPLICATES=true I={input} O={output} \
                        M=output/logs/{wildcards.sample}_marked_dup_metrics.txt 2> {log}")

                if config["keep_unfiltered_bam"] == "FALSE":
                    shell("rm -f {input} {input}.bai")
                ## Rule for calculating insert size for PE sequencing
                shell("picard CollectInsertSizeMetrics \
                    I=output/bam/{wildcards.sample}.sorted.rmdup.bam \
                    O=output/logs/{wildcards.sample}_insert_size_metrics.txt \
                    H=output/logs/{wildcards.sample}_insert_size_histogram.pdf \
                    M=0.05")
                if config["refflat"] != "FALSE":
                    shell("picard CollectRnaSeqMetrics \
                    I=output/bam/{wildcards.sample}.sorted.rmdup.bam \
                    O=output/logs/{wildcards.sample}.RNA_Metrics \
                    REF_FLAT={params.refflat} \
                    STRAND=SECOND_READ_TRANSCRIPTION_STRAND")
            else:
                if config["use_UMI"] == "TRUE":
                    shell("umi_tools dedup --stdin={input} --log={log} > {output}")
                else:
                    shell("cp {input} {output}")
                if config["keep_unfiltered_bam"] == "FALSE":
                    shell("rm -f {input} {input}.bai")
                if config["refflat"] != "FALSE":
                    shell("picard CollectRnaSeqMetrics \
                    I=output/bam/{wildcards.sample}.sorted.bam \
                    O=output/logs/{wildcards.sample}.RNA_Metrics \
                    REF_FLAT={params.refflat} \
                    STRAND=FIRST_READ_TRANSCRIPTION_STRAND")