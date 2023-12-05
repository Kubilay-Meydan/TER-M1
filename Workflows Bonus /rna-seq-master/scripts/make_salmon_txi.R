#!/usr/bin/env Rscript


args = commandArgs(trailingOnly = TRUE)
in_salmon_str = args[1]
out_txi = args[2]
in_txdb = args[3]

# Load libraries
library(paste(in_txdb), character.only=T)
txdb = get(in_txdb)
library(tximport)

# Make salmon file list
in_salmon_vec = unlist(strsplit(in_salmon_str, " "))
names(in_salmon_vec) = substr(gsub("^.*lib", "lib", in_salmon_vec), 1, 6)

# Make gene annotation
k = keys(txdb, keytype = "TXNAME")
tx2gene = AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

# Make txi object
txi = tximport(in_salmon_vec, type = "salmon", tx2gene = tx2gene)

# Save txi object
save(txi, file = out_txi)
