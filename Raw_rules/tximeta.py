## tximeta
rule tximeta:
	input:
	    os.path.join(outputdir, "Rout", "pkginstall_state.txt"),
		expand(os.path.join(outputdir, "salmon", "{sample}", "quant.sf"), sample = samples.names.values.tolist()),
		metatxt = config["metatxt"],
		salmonidx = os.path.join(config["salmonindex"], "versionInfo.json"),
		json = "".join([config["salmonindex"], ".json"]),
		script = "scripts/run_tximeta.R"
	output:
		os.path.join(outputdir, "outputR", "tximeta_se.rds")
	log:
		os.path.join(outputdir, "Rout", "tximeta_se.Rout")
	benchmark:
		os.path.join(outputdir, "benchmarks", "tximeta_se.txt")
	params:
		salmondir = lambda wildcards, input: os.path.dirname(os.path.dirname(input[1])),   ## dirname of second output
		flag = config["annotation"],
		organism = config["organism"],
		Rbin = Rbin
	conda:
		Renv
	shell:
		'''{params.Rbin} CMD BATCH --no-restore --no-save "--args salmondir='{params.salmondir}' json='{input.json}' metafile='{input.metatxt}' outrds='{output}' annotation='{params.flag}' organism='{params.organism}'" {input.script} {log}'''