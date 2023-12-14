rule star_genome:
    input:
        fasta = config['reference']['fasta']['hg38'],
        gtf = config['reference']['gtf']['hg38'],
        sjs =  config['datadirs']['sj_files'] + "/" + "SJ.out.pass1_merged.tab",
        genomedir = config['reference']['stargenomedir']['hg38'],
        queue = rules.merge.output.sjs
    output:
        starindex = config['reference']['stargenomedir']['hg38'] + "/" + "SAindex"
    params:
        overhang = 149
    threads: 12
    resources:
        mem_mb = 40000
    shell: """
        STAR \
        --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {input.genomedir} \
        --outFileNamePrefix {input.genomedir} \
        --genomeFastaFiles {input.fasta} \
        --sjdbGTFfile {input.gtf} \
        --limitSjdbInsertNsj 1637800 \
        --sjdbFileChrStartEnd  {input.sjs} \
        --sjdbOverhang {params.overhang}
          """
