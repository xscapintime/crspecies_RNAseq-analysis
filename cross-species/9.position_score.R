# position score for pathways
# rank based on pathway logFc

rm(list = ls())


library(tidyverse)


## GO
## -------------

load("human_GOallres.Rdata")
load("mouse_GOallres.Rdata")

# go term logFC
read.table("human_bp_fc.txt")
read.table("mouse_bp_fc.txt")

read.table("human_mf_fc.txt")
read.table("mouse_mf_fc.txt")

read.table("human_cc_fc.txt")
read.table("mouse_cc_fc.txt")
