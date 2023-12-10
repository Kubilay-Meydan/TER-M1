rule generate_genome:
	input:
		fasta = get_genome_fasta(config['species']),
		gtf = get_gtf(config['species'])
	output:
		GENOME_DIR + "/SA",
		GENOME_DIR + "/Genome"
	params:
		sjdbOverhang = config['readLen'] - 1
	threads:
		4
	conda:
		"../env/align.yaml"
	shell:
		"""
		STAR \
	    --runThreadN {threads} \
	    --runMode genomeGenerate \
	    --genomeDir {GENOME_DIR} \
	    --genomeFastaFiles {input.fasta} \
	    --sjdbGTFfile {input.gtf} \
	    --sjdbOverhang {params.sjdbOverhang} \
		--genomeSAsparseD 10
		"""
