print("Integration testing snakefile for Post-QC RNA-seq Differential Expression\n")

# Import common packages
import pandas as pd
import re
import numpy as np





rule all:
    input:
        design
        tmm
        ebayes
        dds

rule symlink_rnaseq_de_inputs:
    input:

#include: " <INCLUDE FILE LOCATION (VIA CONFIG PARAM)>"
