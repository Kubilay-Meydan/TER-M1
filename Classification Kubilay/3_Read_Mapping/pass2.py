rule pass2:
   input:
      f1 = config['datadirs']['trim'] + "/" + "{file}_1_val_1.fq.gz",
      f2 = config['datadirs']['trim'] + "/" + "{file}_2_val_2.fq.gz",
      line = rules.star_genome.output.starindex
   output: config['datadirs']['pass2'] + "/" + "{file}_Aligned.toTranscriptome.out.bam", config['datadirs']['pass2'] + "/" + "{file}_Aligned.sortedByCoord.out.bam"
   params:
      genomedir = config['reference']['stargenomedir']['hg38'],
      prefix =  config['datadirs']['pass2'] + "/" + "{file}_"
   threads: 16
   resources:
      mem_mb= 50000
   shell: """
        STAR  \
        --runThreadN {threads} \
        --genomeDir {params.genomedir} \
        --readFilesIn {input.f1} {input.f2} \
        --readFilesCommand zcat \
        --outFileNamePrefix {params.prefix} \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --quantMode TranscriptomeSAM \
        --outSAMattributes NH HI AS NM MD \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --sjdbScore 1 \
        --outBAMsortingThreadN 5 \
        --limitBAMsortRAM 50000000000
        """