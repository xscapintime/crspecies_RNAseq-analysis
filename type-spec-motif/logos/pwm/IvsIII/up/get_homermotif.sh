# TR4(NR),DR1/Hela-TR4-ChIP-Seq(GSE24685)/Homer          1e-18  NA      1e-13  NA
ln -s ../../../../homer/typeI_up_MotifOutput/knownResults/known1.motif TR4.motif

# PPARa(NR),DR1/Liver-Ppara-ChIP-Seq(GSE47954)/Homer     1e-04  NA      1e-05  NA
ln -s ../../../../homer/typeI_up_MotifOutput/knownResults/known21.motif PPARa.motif

# VDR(NR),DR3/GM10855-VDR+vitD-ChIP-Seq(GSE22484)/Homer  0.001  NA      1e-16  NA
ln -s ../../../../homer/typeI_up_MotifOutput/knownResults/known28.motif VDR.motif

# E2F4(E2F)/K562-E2F4-ChIP-Seq(GSE31477)/Homer           0.001  NA      0.001  NA
ln -s ../../../../homer/typeI_up_MotifOutput/knownResults/known33.motif E2F4.motif

# RXR(NR),DR1/3T3L1-RXR-ChIP-Seq(GSE13511)/Homer         NA     1e-190  0.001  NA
wget -c http://homer.ucsd.edu/homer/motif/HomerMotifDB/homerResults/motif251.motif -O RXR.motif

# POL013.1_MED-1/Jaspar                                  NA     1e-63   NA     1e-43
wget -c http://jaspar.genereg.net/api/v1/matrix/POL013.1.pfm -O MED-1.pfm
