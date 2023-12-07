rule DRIMSeq:
	input:
	    os.path.join(outputdir, "Rout", "pkginstall_state.txt"),
		rds = os.path.join(outputdir, "outputR", "edgeR_dge.rds"),
		script = "scripts/run_render.R",
		template = "scripts/DRIMSeq_dtu.Rmd"
	output:
		html = os.path.join(outputdir, "outputR", "DRIMSeq_dtu.html"),
		rds = os.path.join(outputdir, "outputR", "DRIMSeq_dtu.rds")
	params:
		directory = lambda wildcards, input: os.path.dirname(input['rds']),   ## dirname of rds input
		organism = config["organism"],
		ncores = config["ncores"],
		design = config["design"].replace(" ", "") if config["design"] is not None else "",
		contrast = config["contrast"].replace(" ", "") if config["contrast"] is not None else "",
		Rbin = Rbin
	log:
		os.path.join(outputdir, "Rout", "run_dtu_drimseq.Rout")
	benchmark:
		os.path.join(outputdir, "benchmarks", "run_dtu_drimseq.txt")
	conda:
		Renv
	threads:
		config["ncores"]
	shell:
		'''{params.Rbin} CMD BATCH --no-restore --no-save "--args se='{input.rds}' design='{params.design}' contrast='{params.contrast}' ncores='{params.ncores}' rmdtemplate='{input.template}' outputdir='{params.directory}' outputfile='DRIMSeq_dtu.html'" {input.script} {log}'''
