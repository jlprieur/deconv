#!/bin/csh
################################################################
# Reduction procedure for HR6803
# Version 19-01-96 
################################################################
goto bispec
#goto photnoise 
#Syntax is:  
#runs dec_car3 in_photon_file output_file_exten 
#         reduc_fact,total_nphot,phot_per_frame flat_field
# Example: runs dec_car3 test tt 16,1000,100 ff_oc45 
# or without flat field: runs dec_car3 test tt 16,1000,100 0 
# HD33883_f650 
# date =415.2461-08-93
# integration time =606.033020 (sec) 
# nphot =4425484 (photons)
#
#runs decode_cp40 ADS974_oc45_650_1 _1 8,168,140,12,1000,1000000,0 0 0
#runs decode_cp40 ceres_oc45_650_1 _10 8,168,140,12,1000,10000000,0 0 0
#runs decode_cp40 ADS974_oc45_650_1 _2 8,168,140,12,1000,10000000,0 0 0
#runs decode_cp40 SAO80508_oc45_650_1 _20 8,168,140,12,1000,10000000,0 0 0
goto end
photnoise:
### runs photon_corr 12 _5 0 _5p 50.
# in_ext out_ext uv-mask ffield ir_max,max_nclosure phot_per_frame(if ffield=0)
#runs phot_noise_mask _1 _1p 0 0 12,1000 195.6
runs phot_noise_mask _10 _10p 0 0 12,1000 134.8
runs phot_noise_mask _2 _2p 0 0 12,1000 188.1
runs phot_noise_mask _20 _20p 0 0 12,1000 127.0
# modsq_5 is 54x54 pixels, so I get an image of 64x64:
#runs enlarge 1 modsq_5p.fits 64,64 modsq_5c.fits
goto end
bispec:
##############################################################
inv1:
# 2= simulations
# 12= Radius of uv-coverage
# 0= weights set to unity
# 0.08= test to exit from main loop
# 220= number of bispectral terms
# 0= lower cut to take into account spectral terms
# sigma_max sigma_null
rm -f err_mod.dat sigma.dat error1.dat error2.dat bisp_error1.dat
rm -f bisp_error2.dat inv_bispec.log
#modsq_m01c is a subset of modsq_m01
#1.E-4 for modsq_5c or 0.01 for the modulus:
#runs inv_bispec2 1 12,0,0.08,220,0.,0.7,0.5,1000 modsq_1p bisp1_1p
runs inv_bispec2 1 12,0,0.08,220,0.,1.,0.5,1000 modsq_10p bisp1_10p
runs fft2 -1 ref n imf m10 testi
runs inv_bispec2 1 12,0,0.08,220,0.,1.,0.5,1000 modsq_20p bisp1_20p
runs fft2 -1 ref n imf m20 testi
runs inv_bispec2 1 12,0,0.08,220,0.,1.,0.5,1000 modsq_2p bisp1_2p
runs fft2 -1 ref n imf m2 testi
goto end
##############################################################
diane:
##############################################################
# Creates .SNR .SIG .TFR .TFI
##runs d_snr 3. 2. long_200 mm01.
#
# mm01.SIG is set to 1.  
# mm01.FTB is set to 1.  
cp snr1.mt mm01.SNR
cp sigma.mt mm01.SIG
cp ref.mt mm01.TFR
cp imf.mt mm01.TFI
goto end
## Creates mm01.FTB: (not good, so I take a unity file)
##runs recen_fft modsqc_m04 test 
##runs d_ftmm mm01. 0.002 test 
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
runs d_regul $ALPHAT mm01. mask.mt $ALPHAT 
#      generic name, alphat
runs d_erreur mm01. 1.45
runs d_cgradient mm01. 1.,1.,0.
goto end
runs d_fstep mm01. mm01.PHIN 0.04
#rm pstcopy.tmp
#runs pstcopy0 mm01.PHIN
#lpr pstcopy.tmp

**********************************************************************
* Parameters used for datatest.mt
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
