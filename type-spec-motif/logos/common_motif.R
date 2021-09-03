# find common motif between type 1,2,3
# homer and MEME
# ------------------------------------

rm(list = ls())

library(tidyverse)


### homer result
## homer known
homer_known <- read.csv("../homer/typeI_up_MotifOutput/knownResults.txt", header = T, sep = "\t")


## homer de novo