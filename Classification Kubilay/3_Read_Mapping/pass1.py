rule pass1:
   input:
      f1 = config['datadirs']['trim'] + "/" + "{file}_1_val_1.fq.gz",
      f2 = config['datadirs']['trim'] + "/" + "{file}_2_val_2.fq.gz",
      queue = rules.trim_galore_pe.output.rev_pai
   output: config['datadirs']['bam'] + "/" + "{file}_SJ.out.tab", config['datadirs']['bam'] + "/" + "{file}_Aligned.toTranscriptome.out.bam"
   params:
      genomedir = config['reference']['star_ref'],
      prefix =  config['datadirs']['bam'] + "/" + "{file}_"
   threads: 16
   resources:
      mem_mb= 40000
   shell: """
        STAR  \
        --runThreadN {threads} \
        --genomeDir {params.genomedir} \
        --readFilesIn {input.f1} {input.f2} \
        --readFilesCommand zcat \
        --outFileNamePrefix {params.prefix} \
        --outSAMtype None \
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
        --limitBAMsortRAM 50000000000
        """