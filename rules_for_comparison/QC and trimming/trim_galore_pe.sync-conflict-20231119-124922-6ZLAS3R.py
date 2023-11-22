rule trim_galore_pe:
   input:
      f1 = config['datadirs']['fastq'] + "/" + "{file}_1.fq.gz",
      f2 = config['datadirs']['fastq'] + "/" + "{file}_2.fq.gz"
   output:
      fwd_pai = config['datadirs']['trim'] + "/" + "{file}_1_val_1.fq.gz",
      rev_pai = config['datadirs']['trim'] + "/" + "{file}_2_val_2.fq.gz",
   params:
      extra = " -j 8 --illumina -q 20 --phred33 --length 20",  
      prefix =  config['datadirs']['trim'],
   resources:
      mem_mb= 20000
   shell:
        """ 
        trim_galore \
        {params.extra} \
        --paired {input.f1} {input.f2} \
        -o {params.prefix} 
         """   
         