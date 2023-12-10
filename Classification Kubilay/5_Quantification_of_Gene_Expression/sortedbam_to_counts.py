rule sortedbam_to_counts:
    input:
        sorted_bam = "output/bam/{sample}.sorted.bam" if config["type"] == "single" else "output/bam/{sample}.sorted.rmdup.bam"
    output:
        gene_counts = "output/counts/{sample}.gene.counts.txt",
        exon_counts = "output/counts/{sample}.exon.counts.txt"
    params:
        gtf = config["gtf"]
    log:
        "output/logs/{sample}.feature_counts.log"
    run:
        if config["count_scheme"] == "fraction":
            if config["type"] == "paired":
                shell("featureCounts -p -O --fraction  -t gene -a {params.gtf} \
                    -o output/counts/{wildcards.sample}.gene.counts.txt \
                    {input.sorted_bam} 2> output/logs/{wildcards.sample}.gene.feature_counts.log")
                shell("featureCounts -p -O --fraction  -t exon -a {params.gtf} \
                    -o output/counts/{wildcards.sample}.exon.counts.txt \
                    {input.sorted_bam} 2> output/logs/{wildcards.sample}.exon.feature_counts.log")
            if config["type"] == "single":
                shell("featureCounts -O --fraction -t gene \
                    -a {params.gtf} -o output/counts/{wildcards.sample}.gene.counts.txt \
                    {input.sorted_bam} 2> output/logs/{wildcards.sample}.gene.feature_counts.log")
                shell("featureCounts -O --fraction -t exon \
                    -a {params.gtf} -o output/counts/{wildcards.sample}.exon.counts.txt \
                    {input.sorted_bam} 2> output/logs/{wildcards.sample}.exon.feature_counts.log")
        elif config["count_scheme"] == "all_reads":
            if config["type"] == "paired":
                shell("featureCounts -p -O -t gene -a {params.gtf} \
                    -o output/counts/{wildcards.sample}.gene.counts.txt \
                    {input.sorted_bam} 2> output/logs/{wildcards.sample}.gene.feature_counts.log")
                shell("featureCounts -p -O -t exon -a {params.gtf} \
                    -o output/counts/{wildcards.sample}.exon.counts.txt \
                    {input.sorted_bam} 2> output/logs/{wildcards.sample}.exon.feature_counts.log")
            if config["type"] == "single":
                shell("featureCounts -O -t gene -a {params.gtf} \
                    -o output/counts/{wildcards.sample}.gene.counts.txt \
                    {input.sorted_bam} 2> output/logs/{wildcards.sample}.gene.feature_counts.log")
                shell("featureCounts -O -t exon -a {params.gtf} \
                    -o output/counts/{wildcards.sample}.exon.counts.txt \
                    {input.sorted_bam} 2> output/logs/{wildcards.sample}.exon.feature_counts.log")
        elif config["count_scheme"] == "unique_reads":
            if config["type"] == "paired":
                shell("featureCounts -p -t gene -a {params.gtf} \
                    -o output/counts/{wildcards.sample}.gene.counts.txt \
                    {input.sorted_bam} 2> output/logs/{wildcards.sample}.gene.feature_counts.log")
                shell("featureCounts -p -t exon -a {params.gtf} \
                    -o output/counts/{wildcards.sample}.exon.counts.txt \
                    {input.sorted_bam} 2> output/logs/{wildcards.sample}.exon.feature_counts.log")
            if config["type"] == "single":
                shell("featureCounts -t gene -a {params.gtf} \
                    -o output/counts/{wildcards.sample}.gene.counts.txt \
                    {input.sorted_bam} 2> output/logs/{wildcards.sample}.gene.feature_counts.log")
                shell("featureCounts -t exon -a {params.gtf} \
                    -o output/counts/{wildcards.sample}.exon.counts.txt \
                    {input.sorted_bam} 2> output/logs/{wildcards.sample}.exon.feature_counts.log")