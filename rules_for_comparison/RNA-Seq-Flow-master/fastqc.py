rule fastqc:
   input:
      f1 = config['datadirs']['fastq'] + "/" + "{file}_{read}.fq.gz"
   output: config['datadirs']['qc'] + "/" + "{file}_{read}_fastqc.html", config['datadirs']['qc'] + "/" + "{file}_{read}_fastqc.zip"
   params:
      prefix =  config['datadirs']['qc'],
   resources:
      mem_mb= 10000
   shell:
        """
        fastqc  --thread 8 --outdir {params.prefix} --nogroup {input.f1}
        """