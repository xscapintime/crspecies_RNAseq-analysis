install.packages("devtools")
unlink("D:/Program Files/Microsoft/R Open/R-4.0.2/library/00LOCK", recursive = TRUE)

# if in WSL
# forget about it
unlink("/home/shakesbeer/R/x86_64-pc-linux-gnu-library/4.0/00LOCK", recursive = TRUE)


devtools::install_local("D:/Program Files/Microsoft/R Open/R-4.0.2/library/CytoTRACE_0.3.3.tar.gz")




#remove.packages("fs")
#BiocManager::install("sva")

library(CytoTRACE)

Sys.setenv(RETICULATE_PYTHON="C:/Users/Administrator/miniconda3/python.exe") 

## test
results <- CytoTRACE(marrow_10x_expr)

