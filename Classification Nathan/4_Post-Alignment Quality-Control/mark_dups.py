rule mark_dups:
    input:
       bam = config['datadirs']['pass2'] + "/" + "{file}_Aligned.sortedByCoord.out.bam"
    output:
       dbam = config['datadirs']['dedup'] + "/" + "{file}_Aligned.sortedByCoord.out.md.bam",
       metric = config['datadirs']['dedup'] + "/" + "{file}_Aligned.sortedByCoord.out.metrics.txt"
    params:
       picard = "java -jar $EBROOTPICARD/picard.jar"
    resources:
       mem_mb = 10000
    shell: """
         module load picard
        {params.picard} MarkDuplicates  INPUT={input.bam} OUTPUT={output.dbam} METRICS_FILE={output.metric} ASSUME_SORT_ORDER=coordinate OPTICAL_DUPLICATE_PIXEL_DISTANCE=100
          """