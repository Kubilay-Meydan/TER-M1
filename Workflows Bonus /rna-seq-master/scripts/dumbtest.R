#!/usr/bin/env Rscript

############
###      ###
############

# For unit testing


# Command line arguements
args = commandArgs(trailingOnly = TRUE)
out= args[1]

# Load required packages
library(tidyverse)

test = data.frame(top=c(1,2,3),
                  bottom=c('a','b','c'))

test2 = as_tibble(test)

write_tsv(test2, file = out)
