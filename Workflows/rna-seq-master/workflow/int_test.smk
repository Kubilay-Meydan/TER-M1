#########1#########2#########3#########4#########5#########6#########7#########8
###                                                                          ###
###               Integration Testing Snakefile for RNA-seq                  ###
###                                                                          ###
#########1#########2#########3#########4#########5#########6#########7#########8

##################################
###   Load Required Packages   ###
##################################

import numpy as np
import os
import pandas as pd
import re

###########################
###   Variable Naming   ###
###########################

# Names directly from configuration YAML
threads = config['threads']

# Names build from configuration parameter base
rna_script_dir = config['rna_repo'] + "/scripts"

rule dumbtest:
    output: "/tmp/test2.tsv",
    params:
        script = rna_script_dir + "/dumbtest.R"
    shell:
        """
        Rscript {params.script} \
        {output}
        """
