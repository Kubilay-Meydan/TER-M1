rule edgeR:
	input:
		os.path.join(outputdir, "Rout", "pkginstall_state.txt"),
		rds = os.path.join(outputdir, "outputR", "tximeta_se.rds"),
		script = "scripts/run_render.R",
		template = "scripts/edgeR_dge.Rmd"
	output:
		html = os.path.join(outputdir, "outputR", "edgeR_dge.html"),
		rds = os.path.join(outputdir, "outputR", "edgeR_dge.rds")
	params:
		directory = lambda wildcards, input: os.path.dirname(input['rds']),   ## dirname of rds input
		organism = config["organism"],
		design = config["design"].replace(" ", "") if config["design"] is not None else "",
		contrast = config["contrast"].replace(" ", "") if config["contrast"] is not None else "",
		genesets = geneset_param,
		Rbin = Rbin
	log:
		os.path.join(outputdir, "Rout", "run_dge_edgeR.Rout")
	benchmark:
		os.path.join(outputdir, "benchmarks", "run_dge_edgeR.txt")
	conda:
		Renv
	shell:
		'''{params.Rbin} CMD BATCH --no-restore --no-save "--args se='{input.rds}' organism='{params.organism}' design='{params.design}' contrast='{params.contrast}' {params.genesets} rmdtemplate='{input.template}' outputdir='{params.directory}' outputfile='edgeR_dge.html'" {input.script} {log}'''
