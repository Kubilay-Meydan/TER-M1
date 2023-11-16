rule fastqctrim:
   input:
      f1 = config['datadirs']['trim'] + "/" + "{file}_{read}_val_{read}.fq.gz",
      trimmedFiles = rules.trim_galore_pe.output.rev_pai 
   output: config['datadirs']['fatsqctrim'] + "/" +"{file}_{read}_val_{read}_fastqc.html"
   params:
      prefix = config['datadirs']['fatsqctrim'],
   resources:
      mem_mb= 10000
   shell:
        """
       fastqc  --thread 8 --outdir {params.prefix} --nogroup {input.f1}
       """         

