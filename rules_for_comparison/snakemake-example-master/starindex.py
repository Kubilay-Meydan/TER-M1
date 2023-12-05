rule starindex:
	input: ref=FASTAREF, starref=STARREFDIR, 
	output: CHRNAME
	shell: "{STAR} --limitGenomeGenerateRAM 54760833024 --runMode genomeGenerate --genomeDir {input.starref} --genomeFastaFiles {input.ref}"
