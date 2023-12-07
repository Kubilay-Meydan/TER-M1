rule rsem_norm:
   input:
      bam = config['datadirs']['pass2'] + "/" + "{file}_Aligned.toTranscriptome.out.bam"
   output:
      genes = config['datadirs']['quant'] + "/" + "{file}.genes.results"
   params:
      calcexp = config['tools']['rsem']['calcexp'],
      genomedir = config['reference']['rsemgenomedir']['hg38'],
      prefix = config['datadirs']['quant'] + "/" + "{file}"
   threads: 16
   resources:
      mem_mb = 15000
   shell:"""
        {params.calcexp} \
        --paired-end \
        --no-bam-output \
        --quiet \
        --no-qualities \
        -p {threads} \
        --forward-prob 1.0 \
        --seed-length 25 \
        --fragment-length-mean -1.0 \
        --bam {input.bam} {params.genomedir} {params.prefix}
        """