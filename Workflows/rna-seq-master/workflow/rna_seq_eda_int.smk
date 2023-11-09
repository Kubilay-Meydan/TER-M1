#print("Integration testing snakefile for bulk RNA-seq\n")

# Import common packages
import pandas as pd
import re
import numpy as np

datadir = config["datadir"]
inputs=datadir + "/inputs"
analysis = datadir + "/analysis"
salmon = analysis + "/salmon"
results = datadir + "/results"
factor_str= config["factor_str"]
rna_container = config["rna_container"]
logdir = config["datadir"] + "/logs"

rna_repo = config["rna_repo"]
rna_scriptdir = rna_repo + "/scripts"
library_tsv=inputs + "/libraries.tsv"

rna_libraries = pd.read_table(inputs + "/libraries.tsv")
rna_libraries["path"]= inputs + "/" + rna_libraries["basename"]

# Needs full path to work (no tilda)
readable = []
for x in rna_libraries.path:
    readable.append(os.access(x, os.R_OK))
rna_libraries['readable']=readable

rna_libraries = rna_libraries[rna_libraries.readable == True]

rna_library_indict = rna_libraries["library"].tolist()
rna_file_indict = rna_libraries["path"].tolist()
rna_lib_dict = dict(zip(rna_library_indict, rna_file_indict))

BULK_RNA_LIBS = list(rna_lib_dict.keys())

rule all:
    input:
        expand(salmon + "/{library}.quant.sf", library = BULK_RNA_LIBS),
        expand(analysis + "/{experiment}_txi.rdata", experiment = "all"),
        results + "/figures/all_pca.pdf",
        analysis + "/all_eda.rdata",

rule symlink_salmon:
    container: rna_container,
    input: lambda wildcards: rna_lib_dict[wildcards.library],
    log: logdir + "/{library}_symlink_salmon.log",
    output: salmon + "/{library}.quant.sf",
    params:
        script = rna_scriptdir + "/symlink_salmon.sh"
    shell:
        """
        {params.script} {input} {output} &> {log}
        """

include: rna_repo + "/workflow/rna_seq_eda.smk"
