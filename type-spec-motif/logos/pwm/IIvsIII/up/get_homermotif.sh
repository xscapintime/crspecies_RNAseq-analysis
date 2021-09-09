### no overlap in meme results#

# PU.1:IRF8(ETS:IRF)/pDC-Irf8-ChIP-Seq(GSE66899)/Homer  1e-06  NA     1e-06  NA
ln -s ../../../../homer/typeII_up_MotifOutput/knownResults/known6.motif PU.1.motif

# NFkB-p65(RHD)/GM12787-p65-ChIP-Seq(GSE19485)/Homer    1e-05  NA     1e-05  NA
ln -s ../../../../homer/typeII_up_MotifOutput/knownResults/known9.motif NFkB-p65.motif

# HIF-1b(HLH)/T47D-HIF1b-ChIP-Seq(GSE59937)/Homer       1e-04  NA     NA     1e-22
ln -s ../../../../homer/typeII_up_MotifOutput/knownResults/known13.motif HIF-1b.motif

# Zfp281(Zf)/ES-Zfp281-ChIP-Seq(GSE81042)/Homer         1e-04  NA     0.001  NA
ln -s ../../../../homer/typeII_up_MotifOutput/knownResults/known15.motif Zfp281.motif

# E2F3(E2F)/MEF-E2F3-ChIP-Seq(GSE71376)/Homer           0.001  NA     1e-08  1e-54
ln -s ../../../../homer/typeII_up_MotifOutput/knownResults/known19.motif E2F3.motif

# HINFP(Zf)/K562-HINFP.eGFP-ChIP-Seq(Encode)/Homer      0.001  NA     NA     1e-63
ln -s ../../../../homer/typeII_up_MotifOutput/knownResults/known21.motif HINFP.motif

# SF1(NR)/H295R-Nr5a1-ChIP-Seq(GSE44220)/Homer          NA     1e-82  1e-05  NA
wget -c http://homer.ucsd.edu/homer/motif/HomerMotifDB/homerResults/motif253.motif -O SF1.motif

# ZNF652/HepG2-ZNF652.Flag-ChIP-Seq(Encode)/Homer       NA     1e-30  0.001  NA
wget -c http://homer.ucsd.edu/homer/motif/HomerMotifDB/homerResults/motif429.motif -O ZNF652.motif

# POL008.1_DCE_S_I/Jaspar                               NA     1e-17  NA     1e-43
wget -c http://jaspar.genereg.net/api/v1/matrix/POL008.1.pfm -O DCE_S_I.pfm

# Sp1(Zf)/Promoter/Homer                                NA     1e-16  NA     1e-57
wget -c http://homer.ucsd.edu/homer/motif/HomerMotifDB/homerResults/motif344.motif -O Sp1.motif