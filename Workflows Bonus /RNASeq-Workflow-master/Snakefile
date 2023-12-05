import pandas as pd

config_file = "config.yaml"

samples = pd.read_table(config["samples_info", index_col = "sample_name"])

rule all:
    input:
         config["sleuth_out"]

include:"rules/kallisto_index.smk"
include:"rules/kallisto_quant.smk"
include:"rule/sleuth.smk"
