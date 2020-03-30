!****************************************************************** 
! Program in Fortran 90
!
!  spdiv => Direct spectral division
!  wiene => Wiener filter 
!  tikho => Tikhonov regularisation
!  maxen => Maximum Entropy method
!
! JLP
! Version 03/04/2001
!****************************************************************** 
program dcv_deconv_2D
 use dcv_2D_data
 use dcv_2D_mod
 implicit none
 real			:: rms_err,ss,L_phi,L_y,snr
 integer		:: nx0,ny0,iter,pntr
 character(len=1)	:: ans
 character(len=20)	:: prefix
 character(len=60)	:: filename
 character(len=80)	:: comments

call jlp_begin
call jlp_inquifmt

print *,' Program to deconvolve 2D signals'
print *,' JLP version 28/03/2001'

!!!! idim = 128*128 = 16384
!  prefix='s3'
  prefix='s5'

if(simulation)then
  print *,' Simulation: *_psf,*_yy0 are needed (prefix: s3, s5, ...)'
else
  print *,' Real data: *_psf,*_raw are needed'
endif
  print *,' Prefix for filenames: '
  read(5,*) prefix 
  print *,'OK: prefix:',prefix

! Read a file first to determine the size:
if(simulation)then
  filename=trim(prefix)//'_yy0'
else
  filename=trim(prefix)//'_raw'
endif
call jlp_vm_readimag(pntr,nx,ny,filename,comments);
!
! nn = number of pixels of an image (used by frprmn)
nn = nx*ny;
if(nn /= idim)then
  print *,' Fatal error: size of images is not equal to idim=',idim
  stop
endif
! To prepare FFT's:
if(nn > idim)then
  print *,' Fatal error: size of images larger than idim=',idim
  stop
endif
! To prepare FFT's:
print *,' fft_setup with ',nx,ny
call fft_setup(nx,ny)

! Then allocate memory space:
allocate(yy0(nn)); allocate(yy(nn)); allocate(hh(nn));
allocate(yyd(nn)); allocate(w0(nn)); allocate(w1(nn)); 
allocate(w_re(nn)); allocate(w_im(nn));
allocate(hh_re(nn)); allocate(hh_im(nn));
allocate(hht_re(nn)); allocate(hht_im(nn));
allocate(yy_re(nn)); allocate(yy_im(nn));

! Read input files:
!----------------------
if(simulation)then
! Original signal:
  filename=trim(prefix)//'_yy0'
  call jlp_readimag(yy0,nx,ny,nx,filename,comments);
else
! Image to deconvolve 
  filename=trim(prefix)//'_raw'
  call jlp_readimag(yy,nx,ny,nx,filename,comments);
endif

!----------------------
! PSF:
filename=trim(prefix)//'_psf'
call jlp_readimag(hh,nx0,ny0,nx,filename,comments);
! Normalize PSF (to avoid numerical pb...)
hh = hh / maxval(hh)
if(nx /= nx0 .or. ny /= ny0)then 
   print *,' Fatal error: wrong size for psf : nx=',nx0,' ny=',ny0; stop;
endif
! Compute transfer function:
  hh_re = hh; hh_im =0.;
  call fft_float(hh_re,hh_im,nx,ny,1)

!----------------------
! Signal to deconvolve:
if(simulation)then
! Either input it from a file
  if(.false.)then
! Original signal:
    filename=trim(prefix)//'_yyb'
    call jlp_readimag(yy,nx0,ny0,nx,filename,comments);
    if(nx /= nx0 .or. ny /= ny0)then 
       print *,' Fatal error: wrong size for yyb : nx=',nx0,' ny=',ny0; stop;
    endif
  else
! Or create it:  
! snr is in dB
! snr=20 or 13dB
! snr=100 or 20dB : simu 3: with disk
! snr=1000 or 30dB: simu 4: with disk
! snr=20 or 13dB: simu 5: without disk
    snr = 13
    call noisy_signal(snr)
    filename =trim(prefix)//'_yyb'
    write(comments,55) trim(prefix)//': convolution of yy0 with psf, SNR=',snr
55  format(A,F6.2,' dB')
    call jlp_writeimag(yy,nx,ny,nx,filename,comments);
  endif
! End simulation case
endif

! Set to zero all arrays:
w0 = 0.; w1 = 0.; w_re = 0.; w_im = 0.; 

!--------------------------------
! Deconvolution
!--------------------------------

write(6,18); 
18 format(/,&
' ==================== Options: ================',/,&
' - spdiv: spectral division',/,&
' - wiene: Wiener filter',/,&
' Regularisation on x_i only: ',/,&
' - tikho: Tikhonov''s regularisation',/,&
' - ggaus: Generalized Gauss''s regularisation (p=1.1)',/,&
' - maxen: regularisation with maximum entropy',/,&
' - sqrtr: convex sqrt(s2+x2) regularisation',/,&
' Regularisation on (x_i+1 - x_i): ',/,&
' - gmark: Gauss-Markov''s regularisation (p=2)',/,&
' ======== Enter the option you want: ===========')
read(5,*) option
print *,' OK: option=',option
L_phi = 0.; 

positivity = .false.
if(option(1:5) .ne. 'spdiv' .and. option(1:5) .ne. 'wiene'&
   .and. option(1:5) .ne. 'maxen')then
   print *,' Positivity constraint? [y/n]'; read(5,*)ans
   if((ans(1:1) .eq. 'y') .or. (ans(1:1) .eq. 'Y')) positivity=.true.
endif
if(option(1:5) .ne. 'spdiv' .and. option(1:5) .ne. 'wiene')then
   print *,' Enter alpha value: '; read(5,*)alpha
endif

select case(option)

!--------------------------------
! Deconvolution by spectral division:
case('spdiv')
   call dcv_spdiv
   write(comments,10)
10 format('Spectral division:')

!--------------------------------
! Deconvolution by Wiener filter:
case('wiene')
! alpha is the ratio variance of noise / variance noise of signal 
   if(simulation)then
     alpha=10**(-snr/10.)
     print *,'Theoretical value is (sigma_noise/sigma_signal)**2 :',alpha
   endif
   print *,' Enter alpha value: '; read(5,*)alpha
   call dcv_wiener
   write(comments,11) alpha
11 format('Wiener: alpha=',F4.1)

!--------------------------------
! Deconvolution with Tikhonov's regularisation 
! Generalized Gauss with p=2
case('tikho')
   call dcv_tikhonov(iter)
   write(comments,12) alpha,iter
12 format('Tikhonov: alpha=',F4.1,' iter=',I3)
   L_phi = phi_tikho(yyd)

!--------------------------------
! Deconvolution with generalized Gauss' regularisation 
! with p=1.1
case('ggaus')
   call dcv_ggauss(iter)
   write(comments,13) alpha,iter
13 format('Gen. Gauss: alpha=',F4.1,' it=',I3)
   L_phi = phi_ggauss(yyd)

!--------------------------------
! Deconvolution with convex sqrt(s2+x2) regularisation 
case('sqrtr')
   print *,' sqrt(s^2+x^2)/Enter s value: '; read(5,*)ss
   call dcv_sqrtrf(ss,iter)
   write(comments,14) ss,alpha,iter
14 format('sqrt: s=',F5.3,' alpha=',F4.1,' it=',I3)
   L_phi = phi_sqrtrf(yyd)

!--------------------------------
! Deconvolution with Maximum Entropy method 
case('maxen')
   call dcv_mem(iter)
   write(comments,15) alpha,iter
15 format('MEM: alpha=',F4.1,' iter=',I3)
   L_phi = phi_mem(yyd)

!--------------------------------
! Deconvolution with Gauss Markov regularisation (on x_i+1 - x_i) 
case('gmark')
   call dcv_gmark(iter)
   write(comments,16) alpha,iter
16 format('Gauss-Markov alpha=',F4.1,' it=',I3)
   L_phi = phi_gmark(yyd)

!---------------------------------
case default
   print *,' Fatal: invalid option'
   stop
end select

if(simulation)then
   rms_err=sqrt(sum((yyd-yy0)**2))/float(nn)
   write(comments,22) trim(comments),rms_err
22 format(a,' rms=',e10.3)
endif

!----------------- L Curve:
! L_phi = phi(x)
! L_y = || y - H x ||^2
    if(L_phi /= 0.)then
!      alpha = 0.;
      L_y = func(yyd)
      write(6,33) L_phi,L_y 
33    format('L curve: L_phi=',G10.4,' L_y=',G10.4)
    endif

filename =trim(prefix)//'_yyd'
call recent_fft(yyd,yyd,nx,ny,nx)
call jlp_writeimag(yyd,nx,ny,nx,filename,comments);


call jlp_end

end program dcv_deconv_2D 
!----------------------------------------------------------------
! Generates a noisy signal with a given SNR
! snr: in dB 
!----------------------------------------------------------------
subroutine noisy_signal(snr)
 use dcv_2D_data
 use dcv_2D_mod
 implicit none
real			:: snr,noise_level
!
! dcv_conv_hh returns: yy = conv(hh,yy0)
  call dcv_conv_hh(yy0,hh_re,hh_im,yy)
  call recent_fft(yy,yy,nx,ny,nx)
! Generate noisy signal with SNR in dB:
! Multiply by 12, since variance of uniform random function is 1/12.
  noise_level = 12. * (sum(yy**2) * 10**(-snr/10.)) / float(nn)
  noise_level = sqrt(noise_level)
  write(6,'("Simulated signal with SNR=",f8.2,"dB => noise level:",f8.3)') snr,noise_level
  call random_seed; call random_number(w0);
! Center noise:
  w0 = w0 - 0.5;
!
  yy = yy + w0 * noise_level;
return
end
