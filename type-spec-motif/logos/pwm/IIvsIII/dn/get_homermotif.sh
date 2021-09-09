## no overlap in meme results

# Oct4:Sox17(POU,Homeobox,HMG)/F9-Sox17-ChIP-Seq(GSE44553)/Homer    1e-07  NA      1e-05  NA
ln -s ../../../../homer/typeII_dn_MotifOutput/knownResults/known5.motif Oct4.motif

# Hoxa10(Homeobox)/ChickenMSG-Hoxa10.Flag-ChIP-Seq(GSE86088)/Homer  1e-04  NA      1e-09  NA
ln -s ../../../../homer/typeII_dn_MotifOutput/knownResults/known15.motif Hoxa10.motif

# PRDM1/MA0508.2/Jaspar                                             NA     1e-136  NA     1e-74
wget -c http://jaspar.genereg.net/api/v1/matrix/MA0508.2.pfm -O PRDM1.motif

# MITF(bHLH)/MastCells-MITF-ChIP-Seq(GSE48085)/Homer                NA     1e-85   NA     1e-53
wget -c http://homer.ucsd.edu/homer/motif/HomerMotifDB/homerResults/motif216.motif -O MITF.motif
