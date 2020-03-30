#!/bin/csh
################################################################
# Reduction procedure for HR6803
# Deconvolution
# Version 19-01-96 
################################################################
#goto inv1
goto diane
##############################################################
inv1:
# 2= simulations
# 12= Radius of uv-coverage
# 0= weights set to unity
# 0.008= test to exit from main loop
# 220= number of bispectral terms
# 0= lower cut to take into account spectral terms
#rm err_mod.dat error1.dat error2.dat bisp_error1.dat bisp_error2.dat
rm sigma.dat inv_bispec2.log
runs inv_bispec2 1 12,10,0.008,120,0.,1.0,0.5,2000 modsq_2345p bisp1_2345p
runs fft2 -1 ref n imf V586 testi
goto end
##############################################################
diane:
##############################################################
# Creates .SNR .SIG .TFR .TFI
##runs d_snr 3. 2. long_200 mm01.
#
# mm01.SIG is set to 1.  
# mm01.FTB is set to 1.  
cp snr1.fits mm01.SNR
cp sigma.fits mm01.SIG
cp ref.fits mm01.TFR
cp imf.fits mm01.TFI
goto end
## *.FTB: modulus of the limited transfer function.
## Creates mm01.FTB from modsq_*:
## (zero freq. is at the central pixel)
runs d_ftmm mm01. 0.002 modsq_112p
goto end
# Creates .HR .FLI .PHIT .PCR .PCI
# alphaT, generic name, Hr radius
begin:
setenv ALPHAT 1.42
runs d_lissage $ALPHAT mm01. 8.
#
# Plot first approximation of the solution
#rm pstcopy.tmp
#runs pstcopy0 mm01.PHIT
#lpr pstcopy.tmp
#       isophote or mask_file, alphat
#
# Add support constraint:
#runs d_regul 4. mm01. 3.E-4 1.4  
runs d_regul $ALPHAT mm01. mask.fits $ALPHAT 
#      generic name, alphat
runs d_erreur mm01. 1.45
runs d_cgradient mm01. 1.,1.,0.
goto end
runs d_fstep mm01. mm01.PHIN 0.04
#rm pstcopy.tmp
#runs pstcopy0 mm01.PHIN
#lpr pstcopy.tmp

**********************************************************************
* Parameters used for datatest.fits
**********************************************************************

 D_SNR:  snr=10. reliability=2
 [[ D_REPONS:   psf : sigma=3  ---> rayon=12. ]]
 D_FTM: lim = lower threshold of the PSF in Fourier domain 
 D_LISSAGE: gain=2, alphaT=1.5 psf radius (in direct space)=35.
            ---> eta0 (=1.76)
 D_REGUL: isophote threshold=112., alphaprime=4.
 D_ERREUR: ---> epsilon (=586.,586.,0.)
 D_CGRADIENT: epsilon ---> *.PHIN (final image)

**********************************************************************
* Input/Output
**********************************************************************

D_SNR : 
        IN: input image
        OUT: *.SNR (SNR of the input image, Fourier domain) 
             *.SIG (noise of FT), *.TFR (FT real of the input)
-             *.TFI (FT imagin. of the input)

D_REPONS : 
        OUT: *.PSF (Point spread function, direct space)

D_FTM : 
        IN: *.PSF
        OUT: *.FTB (modulus of the limited transfer function)

D_LISSAGE: 
        IN: *.FTB, *.SNR, *.TFR, *.TFI
        OUT: *.HR (synthetic aperture), *.FLI (smoothing filter)
             *.PHIT  (first approx. of the object),
             *.PCR and *.PCI (real and imag. parts of the FT of *.PHIT)

D_REGUL:
        IN: *.PHIT, *.SNR, *.HR
        OUT: *.V (Domaine of the object, direct space), *.G (regularizing
             function), *.KR (Stabilizing function)
D_ERREUR:
        IN: *.G, *.FLI, *.FTB, *.SNR, *.SIG, *.KR

D_CGRADIENT:
        IN: *.PHIT, *.PCR, *.PCI, *.V, *.G
        OUT: *.PHIN (final estimation of the object)

end:
