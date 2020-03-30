!****************************************************************** 
! Program in Fortran 90
! To emulate Herve Carfantan's simulations PI1 and PI2 
!
! JLP
! Version 19/03/2001
!****************************************************************** 
program dcv_test1 
 implicit none
 integer,parameter		:: nn1=256
 logical,parameter		:: display_all=.false.
 logical        		:: simu1
 integer			:: nn,i,hh_nx
 real, dimension(nn1)		:: xx,yy0,yy0_pw,hh,yy,yyb,hh_pw,w0,w1
 real, dimension(nn1)		:: yy_noise 
! snr in dB (20 db = 100.)
 real,parameter			:: snr=20.
 real				:: noise_level
 character(len=1)		:: ans 
 character(len=32)		:: plotdev
 character(len=40)		:: xlabel,ylabel,title 
 character(len=20)		:: extension 

call jlp_begin
call jlp_inquifmt

print *,' Program to generate a signal for simulations'
print *,' JLP version 22/03/2001'
print *,' Menu: '
print *,' 1. Simulation PI1'
print *,' 2. Simulation PI2'
print *,' Enter your choice:'
read(5,10) ans
10 format(a1)

if(ans.eq.'1')then
  simu1=.true.
  print *,' Simulation PI1'
  nn=100
else
  simu1=.false.
  print *,' Simulation PI2'
  nn=250
endif

xlabel=" "; ylabel=" "; plotdev="&xterm";
! Set to zero all arrays:
yy0 = 0.; hh = 0.; yy = 0.; yyb = 0.;

! To prepare FFT's:
call fft_setup(nn1,1)

! X axis:
xx(1:nn1) = (/ (i,i=1,nn1) /);

! Generate Herve's signal
if(simu1)then
  call herve_signal1(yy0,nn)
else
  call herve_signal2(yy0,nn,nn1)
endif

! Display original signal:
if(display_all)then
  write(6,*) "Display original signal:"
  title="Original signal"; 
  call display1(xx,yy0,1,nn,xlabel,ylabel,title,plotdev)
endif

! Compute power spectrum of original signal:
  yy0_pw=yy0; w0=0.;
  call fft_float(yy0_pw,w0,nn1,1,1);
print *,' Sum(signal)',sum(yy0)
print *,' Central value of Real part of FT:',yy0_pw(1)
print *,' Ratio:',sum(yy0)/yy0_pw(1)
  yy0_pw = yy0_pw*yy0_pw + w0*w0;
  call recent_fft_1d_x_float(yy0_pw,yy0_pw,nn1,1,nn1);

! Display power spectrum of original signal:
if(display_all)then
  write(6,*) "Display power spectrum of original signal:"
  title="Power spectrum of original signal" 
  w0 = (/ (float(i)/float(nn1),i=-nn1/2,nn1/2-1) /)
  w1 = log10(yy0_pw)
  call display1(w0,w1,1,nn1,xlabel,ylabel,title,plotdev)
endif

! Generate Herve's PSF:
if(simu1)then
  call herve_psf1(hh,nn)
  hh_nx=nn
else
  call herve_psf2(hh,nn,hh_nx)
endif

! Display PSF:
if(display_all)then
  write(6,*) "Display PSF of filter:"
  title="PSF of the filter" 
  call display1(xx,hh,1,hh_nx,xlabel,ylabel,title,plotdev)
endif

! Compute transfer function:
hh_pw = hh; w0 = 0.;
! Arguments: (real,imag,nx,ny,direct=1 or inverse=-1)
call fft_float(hh_pw,w0,nn1,1,1);
hh_pw = hh_pw*hh_pw + w0*w0;
call recent_fft_1d_x_float(hh_pw,hh_pw,nn1,1,nn1);

! Display Transfer function:
if(display_all)then
  write(6,*) "Frequency response of the filter:"
  title="Frequency response of the filter" 
  w0 = (/ (float(i)/float(nn1),i=-nn1/2,nn1/2-1) /)
  w1 = log10(hh_pw)
  call display1(w0,w1,1,nn1,xlabel,ylabel,title,plotdev)
endif

! Generate filtered signal:
call dcv_conv1D(yy0,hh,yy,nn1,"full")

! Display filtered signal:
if(display_all)then
  write(6,*) "Display filtered signal:"
  title="Signal filtered by the PSF" 
  call display1(xx,yy,1,nn,xlabel,ylabel,title,plotdev)
endif

! Compute power spectrum of filtered signal: 
if(display_all)then
  w0 = yy; w1 = 0.;
  call fft_float(w0,w1,nn1,1,1);
  w0 = w0*w0 + w1*w1;
  call recent_fft_1d_x_float(w0,w0,nn1,1,nn1);
  write(6,*) "Display power spectrum of filtered signal:"
  title="Power spectrum of filtered signal" 
  w1 = log10(w0)
  w0 = (/ (float(i)/float(nn1),i=-nn1/2,nn1/2-1) /)
  call display1(w0,w1,1,nn1,xlabel,ylabel,title,plotdev)
endif

! Generate noisy signal:
! Multiply by 12, since variance of uniform random function is 1/12:
noise_level = 12 * (sum(yy**2) * 10**(-snr/10.)) / float(nn)
noise_level = sqrt(noise_level) 
write(6,'("Simulated signal with SNR=",f8.2,"dB => noise level:",f8.3)') snr,noise_level
call random_seed; call random_number(yy_noise);
! Center the noise (since mean=0.5)
yy_noise = yy_noise - 0.5;
yyb = yy + yy_noise * noise_level;

! Display the simulated (noisy and filtered) signal:
if(display_all)then
  write(6,*) "Display simulated signal:"
  title="Filtered signal with noise (SNR=20)" 
  call display1(xx,yyb,1,nn,xlabel,ylabel,title,plotdev)
endif

! Compute power spectrum of filtered signal: 
!if(display_all)then
if(.true.)then
  w0 = yyb; w1 = 0.;
  call fft_float(w0,w1,nn1,1,1);
  w0 = w0*w0 + w1*w1;
  call recent_fft_1d_x_float(w0,w0,nn1,1,nn1);
  write(6,*) "Display power spectrum of noisy filtered signal:"
  title="Pow. sp. of noisy filtered signal" 
  w1 = log10(w0)
  w0 = (/ (float(i)/float(nn1),i=-nn1/2,nn1/2-1) /)
  call display1(w0,w1,1,nn1,xlabel,ylabel,title,plotdev)
endif

if(simu1)then
  extension='simu1D'
else
  extension='simu2D'
endif
call output_fits(yy0,yy,yyb,hh,nn,extension) 

call jlp_end

end program dcv_test1 
!************************************************
! To emulate Herve Carfantan's simulation PI1 
!***********************************************
subroutine herve_signal1(yy0,nn)
 implicit none
 integer                        :: nn,i
 real, dimension(nn)            :: yy0

! Y axis (Cf. Herve Carfantan's simulation):
  yy0(1:10) = 1.;
  yy0(11:26) = (/ (1.+0.25*float(i), i=1,16) /);
  yy0(27:46) = 2.;
  yy0(47:67) = (/ (2.-0.1*float(i-1), i=1,21) /);
  yy0(68:77) = 4.;
  yy0(78:nn) = 0.;

! Debug:
!  write(6,*) (yy0(i),i=1,nn)

return
end
!************************************************
! To emulate Herve Carfantan's simulation PI1 
! PSF of the filter (rectangle) 
!***********************************************
subroutine herve_psf1(hh,nn)
 implicit none
 integer			:: nn
 real, dimension(nn)		:: hh 

! Y axis (Cf. Herve Carfantan's simulation): 
  hh(1:5) = 1.;
  hh(6:nn) = 0.;

return
end
!************************************************
! To emulate Herve Carfantan's simulation PI2 
!***********************************************
subroutine herve_signal2(yy0,nn,nn1)
 implicit none
 integer                        :: nn,nn1
 real, dimension(nn1)           :: yy0
 integer			:: i

print *,'herve_signal2: simulation #2  nn=',nn,' nn1=',nn1
open(unit=1,file="bgg.asc",status="old")

read(1,*,err=99) (yy0(i),i=1,nn)
yy0(nn+1:nn1)=0.
!
!format(2X,E14.7)

! Debug:
! print *,' bgg signal '
! write(6,*) (yy0(i),i=1,nn)

close(1)
return
! Error case:
99 print *,'Fatal error reading bgg.asc'
   write(6,*) (yy0(i),i=1,nn)
close(1)
stop
end
!************************************************
! To emulate Herve Carfantan's simulation PI2 
! PSF of the filter (rectangle) 
!***********************************************
subroutine herve_psf2(hh,nn,hh_nx)
 implicit none
 integer			:: nn,hh_nx
 real, dimension(nn)		:: hh 

print *,' signal #2 '
open(unit=1,file="ricker_ri.asc",status="old")

hh_nx=32
 read(1,*) hh(1:hh_nx)
 hh(hh_nx+1:256)=0.

return
end
!************************************************
! To emulate Herve Carfantan's conv1D
! IN:
!  yy0: signal to filter
!  hh: filter
!  nn: size of the original signal
!  option: "full", "pre ", "post"
! OUT:
!  yy: filtered signal
!***********************************************
subroutine dcv_conv1D(yy0,hh,yy,nn,option)
 implicit none
 integer                        :: nn
 real, dimension(nn)            :: yy0,hh,yy
 character(len=4)               :: option

select case(option)
   case("full")
     call dcv_conv_1D(yy0,hh,yy,nn)
   case("pre ")
   case("post")
   case default
    print *,' dcv_conv1D/Fatal: bad option'
    stop
end select

return
end
!************************************************
! Convolution product
! IN:
!  yy0: signal to filter
!  hh: filter
!  nn: size of the yy0, hh, arrays (here should be equal to nn1=256) 
!  option: "full", "pre ", "post"
! OUT:
!  yy: filtered signal
!***********************************************
subroutine dcv_conv_1D(yy0,hh,yy,nn)
 implicit none
 integer                        :: nn
 real, dimension(nn)            :: yy0,hh,yy
 real, dimension(nn) 		:: yy0_re,yy0_im,hh_re,hh_im,ww
 character(len=4)		:: option

yy0_re = yy0; yy0_im = 0.;
call fft_float(yy0_re,yy0_im,nn,1,1)
hh_re = hh; hh_im = 0.;
call fft_float(hh_re,hh_im,nn,1,1)
yy = yy0_re*hh_re - yy0_im*hh_im;
ww = yy0_re*hh_im + yy0_im*hh_re;
call fft_float(yy,ww,nn,1,-1)

return
end
!***********************************************
!
!***********************************************
subroutine output_fits(yy0,yy,yyb,hh,nn,extension) 
 implicit none
 integer                        :: nn
 real, dimension(nn) 		:: yy0,yy,yyb,hh
 character(len=20)		:: extension 
 character(len=60)		:: filename 
 character(len=80)		:: comments 
 integer(kind=4)                :: nx,ny

nx=nn; ny=1;
filename=trim(extension)//'_yy0'
comments='Original signal'
call jlp_writeimag(yy0,nx,ny,nx,filename,comments);

filename=trim(extension)//'_yy'
comments='Filtered signal'
call jlp_writeimag(yy,nx,ny,nx,filename,comments);

filename=trim(extension)//'_yyb'
comments='Simulated signal'
call jlp_writeimag(yyb,nx,ny,nx,filename,comments);

filename=trim(extension)//'_psf'
comments='PSF'
call jlp_writeimag(hh,nx,ny,nx,filename,comments);

return
end

