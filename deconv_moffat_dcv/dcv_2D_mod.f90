!****************************************************************** 
! Module used by dcv_deconv_2D.f90
!
!  spdiv => Direct spectral division
!  wiene => Wiener filter 
!  tikho => Tikhonov regularisation
!  maxen => Maximum Entropy method
!
! JLP
! Version 28/03/2001
!****************************************************************** 
module dcv_2D_mod
use dcv_2D_data
implicit none
 private

 public		:: dcv_wiener,dcv_spdiv,dcv_mem,dcv_sqrtrf
 public		:: dcv_tikhonov,dcv_ggauss,dcv_gmark
 public		:: dcv_conv_hh,func,dfunc
! Needed for L curve:
 public		:: phi_ggauss,phi_gmark,phi_tikho,phi_sqrtrf,phi_mem
 private	:: dphi_ggauss,dphi_gmark,dphi_tikho,dphi_sqrtrf,dphi_mem
 private 	:: dcv_check_grad

!------------------
! For dcv_2D_mod:
! ss2 is a parameter needed for convex sqrt regularization
 real,private			:: ss2
! used for calibrating constants 
 real,private			:: gscale=1.e6 

!------------------
! For dcv_cgrad_mod:
 integer,private		:: ncom
 real,private,dimension(idim)	:: pcom,xicom
 public				:: frprmn
 private			:: norm_l2,linmin,df1dim,f1dim,mnbrak,dbrent

contains
!*************************************************************
! Deconvolution with Gauss-Markov regularisation 
! Criterium to minimize is:
! || y - H x ||^2 + alpha phi
!
! with: 
!        _
!        \                     p                      p
! phi =  / | x        - x     |  + | x       - x     |
!        -    (i+1,j)    (i,j)        (i,j+1)   (i,j)
!        i
! (Smooth function, which reduces the variation between two successive pixels)
!*************************************************************
subroutine dcv_gmark(iter)
 integer        	:: iter,i,j
 real			:: ftol,fret

! Compute transfer function:
hh_re = hh; hh_im = 0.;
call fft_float(hh_re,hh_im,nx,ny,1);

! Compute FFT function of reversed PSF (i.e., PSF(-x)):
do j=1,ny
  do i=1,nx
   hht_re(i + j*nx) = hh((nx - (i-1)) + (ny - (j-1))*nx)
  end do
end do
hht_im = 0.;
call fft_float(hht_re,hht_im,nx,ny,1);

! Check if gradient is OK
!call dcv_check_grad(yyd,w0,w1,nn)

! Starting point:
    yyd = 0. 
! Tolerance
    ftol=1.e-6
    call frprmn(yyd,nn,ftol,iter,fret)
    write(6,11) iter,fret
11  format(" Number of iterations:",I4,/,&
           " Value of the minimum:",E12.5)
return
end subroutine dcv_gmark
!*************************************************************
! Regularisation function to be minimized for Gauss-Markov 
!        _
!        \                     p                      p
! phi =  / | x        - x     |  + | x       - x     |
!        -    (i+1,j)    (i,j)        (i,j+1)   (i,j)
!        i
! (Smooth function, which reduces the variation between two successive pixels)
!*************************************************************
real function phi_gmark(xx)
 real,dimension(nx,ny)	:: xx
 double precision	:: dx,dy,ssum 
 integer		:: i,j

ssum = 0.
! Central part:
do j=1,ny-1
  do i=1,nx-1
     dx = xx(i+1,j) - xx(i,j)
     dy = xx(i,j+1) - xx(i,j)
     ssum = ssum + dx*dx + dy*dy
  end do
end do
! Wrap the image for the edges
! Last line
do i=1,nx-1
     dx = xx(i+1,ny) - xx(i,ny)
     dy = xx(i,1) - xx(i,ny)
     ssum = ssum + dx*dx + dy*dy
end do
! Last column 
do j=1,ny-1
     dx = xx(1,j) - xx(nx,j)
     dy = xx(nx,j+1) - xx(nx,j)
     ssum = ssum + dx*dx + dy*dy
end do
! Last pixel
     dx = xx(1,ny) - xx(nx,ny)
     dy = xx(nx,1) - xx(nx,ny)
     ssum = ssum + dx*dx + dy*dy

phi_gmark = real(ssum) 

return
end function phi_gmark
!*************************************************************
! Gradient of potential function for Gauss-Markov 
!        _
!        \                     p                      p
! phi =  / | x        - x     |  + | x       - x     |
!        -    (i+1,j)    (i,j)        (i,j+1)   (i,j)
!        i
! (Smooth function, which reduces the variation between two successive pixels)
!
!*************************************************************
subroutine dphi_gmark(xx,dx)
 integer		:: i,j
 real,dimension(nx,ny)	:: xx,dx

! Gradient is linear, so
! First go along the lines
do i=2,nx-1
  do j=1,ny
   dx(i,j) = - 2*(xx(i+1,j)-xx(i,j)) + 2*(xx(i,j)-xx(i-1,j))
  end do
end do
! First and last pixels of each line:
do j=1,ny
   dx(1,j) = - 2*(xx(2,j)-xx(1,j)) + 2*(xx(1,j)-xx(nx,j))
   dx(nx,j) = - 2*(xx(1,j)-xx(nx,j)) + 2*(xx(nx,j)-xx(nx-1,j))
end do

! Then go along the columns:
do i=1,nx
  do j=2,ny-1
   dx(i,j) = dx(i,j) - 2*(xx(i,j+1)-xx(i,j)) + 2*(xx(i,j)-xx(i,j-1))
  end do
end do
! First and last pixels of each column:
do i=1,nx
   dx(i,1) = dx(i,1) - 2*(xx(i,2)-xx(i,1)) + 2*(xx(i,1)-xx(i,ny))
   dx(i,ny) = dx(i,ny) - 2*(xx(i,1)-xx(i,ny)) + 2*(xx(i,ny)-xx(i,ny-1))
end do

 return
end subroutine dphi_gmark
!*************************************************************
! Deconvolution with MEM 
! Criterium to minimize is:
! || y - H x ||^2 + alpha Sigma_{k=1}^{k=N} x_k^2
!*************************************************************
subroutine dcv_mem(iter)
 integer        	:: iter,i,j
 real			:: ftol,fret

! Compute transfer function:
hh_re = hh; hh_im = 0.;
call fft_float(hh_re,hh_im,nx,ny,1);

! Compute FFT function of reversed PSF (i.e., PSF(-x)):
do i=1,nx
  do j=1,ny
   hht_re(i + j*nx) = hh((nx - (i-1)) + (ny - (j-1))*nx)
  end do
end do
hht_im = 0.;
call fft_float(hht_re,hht_im,nx,ny,1);

! Starting point:
    yyd = 0. 
! Tolerance
    ftol=1.e-6
    call frprmn(yyd,nn,ftol,iter,fret)
    write(6,11) iter,fret
11  format(" Number of iterations:",I4,/,&
           " Value of the minimum:",E12.5)
return
end subroutine dcv_mem
!*************************************************************
! Regularisation function to be minimized for MEM 
! Entropy potential function:
!        _
!        \     
! phi =  /   x log(x )
!        -    i     i
!        i
!*************************************************************
real function phi_mem(xx)
 real,dimension(idim)	:: xx

where(xx > 0.) 
  w1 = xx*log(xx);
elsewhere
  w1 = 1000.*gscale;
endwhere

phi_mem = sum(w1)

return
end function phi_mem
!*************************************************************
! Gradient of potential function for MEM 
!        _
!        \     
! phi =  /   x log(x )
!        -    i     i
! => dphi = sum(1 + log(x)):
!*************************************************************
subroutine dphi_mem(xx,dx)
 real,dimension(idim)	:: xx,dx

where(xx > 0.) 
  dx = 1. + log(xx);
elsewhere
! Next statement is wrong:
!  dx = -1000.*gscale;
! So I leave it:
  dx = -1000.;
endwhere

 return
end subroutine dphi_mem
!*************************************************************
! Deconvolution by Tikhonov's regularisation 
! Criterium to minimize is:
! || y - H x ||^2 + alpha Sigma_{k=1}^{k=N} x_k^2
! (Carfantan: grf)
!*************************************************************
subroutine dcv_tikhonov(iter)
 integer        	:: iter,i,j
 real			:: ftol,fret

! Compute transfer function:
hh_re = hh; hh_im = 0.;
call fft_float(hh_re,hh_im,nx,ny,1);

! Compute FFT function of reversed PSF (i.e., PSF(-x)):
do j=1,ny
  do i=1,nx
   hht_re(i + j*nx) = hh((nx - (i-1)) + (ny - (j-1))*nx)
  end do
end do
hht_im = 0.;
call fft_float(hht_re,hht_im,nx,ny,1);

! Check if gradient is OK
!call dcv_check_grad(yyd,w0,w1,nn)

! Starting point:
    yyd = 0. 
! Tolerance
    ftol=1.e-6
    call frprmn(yyd,nn,ftol,iter,fret)
    write(6,11) iter,fret
11  format(" Number of iterations:",I4,/,&
           " Value of the minimum:",E12.5)
return
end subroutine dcv_tikhonov
!*************************************************************
! Regularisation function to be minimized for Tikhonov 
! Bauman & Sauer's potential function, with p=2
!        _
!        \        p
! phi =  /  | x  |
!        -     i
!        i
!  phi = sum(abs(real(x(:))).^p);
!  dphi = p*sign(x).*(abs(real(x)).^(p-1));
!*************************************************************
real function phi_tikho(xx)
 real,dimension(idim)	:: xx

phi_tikho = sum(xx**2) 

return
end function phi_tikho
!*************************************************************
! Gradient of potential function for Tikhonov's regularisation:
! Bauman & Sauer's potential function, with p=2
!        _
!        \        p
! phi =  /  | x  |
!        -     i
!        i
!  dphi = p*sign(x).*(abs(real(x)).^(p-1));
! here: dphi = 2 x
!
!*************************************************************
subroutine dphi_tikho(xx,dx)
 real,dimension(idim)	:: xx,dx

 dx = 2*xx

 return
end subroutine dphi_tikho
!*************************************************************
! Deconvolution by generalized Gauss' regularisation 
! Criterium to minimize is:
! || y - H x ||^2 + alpha Sigma_{k=1}^{k=N} x_k^p
! with p=1.1
! (Carfantan: ggrf)
!*************************************************************
subroutine dcv_ggauss(iter)
 integer        	:: iter,i,j
 real			:: ftol,fret

! Compute transfer function:
hh_re = hh; hh_im = 0.;
call fft_float(hh_re,hh_im,nx,ny,1);

! Compute FFT function of reversed PSF (i.e., PSF(-x)):
do j=1,ny
  do i=1,nx
   hht_re(i + j*nx) = hh((nx - (i-1)) + (ny - (j-1))*nx)
  end do
end do
hht_im = 0.;
call fft_float(hht_re,hht_im,nx,ny,1);

! Starting point:
    yyd = 0. 
! Tolerance
    ftol=1.e-6
    call frprmn(yyd,nn,ftol,iter,fret)
    write(6,11) iter,fret
11  format(" Number of iterations:",I4,/,&
           " Value of the minimum:",E12.5)
return
end subroutine dcv_ggauss
!*************************************************************
! Regularisation function to be minimized for generalized Gauss
! Bauman & Sauer's potential function, with p=1.1
!        _
!        \        p
! phi =  /  | x  |
!        -     i
!        i
!  phi = sum(abs(real(x(:))).^p);
!  dphi = p*sign(x).*(abs(real(x)).^(p-1));
!*************************************************************
real function phi_ggauss(xx)
 real,dimension(idim)	:: xx

phi_ggauss = sum(abs(xx)**1.1)

return
end function phi_ggauss
!*************************************************************
! Gradient of potential function for generalized Gauss' regularisation:
! Bauman & Sauer's potential function, with p=1.1
!        _
!        \        p
! phi =  /  | x  |
!        -     i
!        i
!  dphi = p*sign(x).*(abs(real(x)).^(p-1));
! here: dphi = 1.1 x**0.1
!
!*************************************************************
subroutine dphi_ggauss(xx,dx)
 real,dimension(idim)	:: xx,dx

where(xx > 0.)
  w1 = xx**0.1
elsewhere
  w1 = -((-xx)**0.1)
end where
 dx = 1.1*w1

 return
end subroutine dphi_ggauss
!*************************************************************
! Deconvolution by generalized Gauss' regularisation 
! Criterium to minimize is:
!             || y - H x ||^2 + alpha phi 
! with:
!        _      __________
!        \     /  2     2 |
! phi =  / \  / ss  +  x
!        -  \/           i
!        i
! (Carfantan: sqrtrf, with ss=0.01 or ss=0.001)
!*************************************************************
subroutine dcv_sqrtrf(ss,iter)
 integer        	:: iter,i,j
 real			:: ss,ftol,fret

! Set ss2 (private parameter) to the value chosen when calling the routine:
ss2 = ss*ss

! Compute transfer function:
hh_re = hh; hh_im = 0.;
call fft_float(hh_re,hh_im,nx,ny,1);

! Compute FFT function of reversed PSF (i.e., PSF(-x)):
do j=1,ny
  do i=1,nx
   hht_re(i + j*nx) = hh((nx - (i-1)) + (ny - (j-1))*nx)
  end do
end do
hht_im = 0.;
call fft_float(hht_re,hht_im,nx,ny,1);

! Starting point:
    yyd = 0. 
! Tolerance
    ftol=1.e-6
    call frprmn(yyd,nn,ftol,iter,fret)
    write(6,11) iter,fret
11  format(" Number of iterations:",I4,/,&
           " Value of the minimum:",E12.5)
return
end subroutine dcv_sqrtrf
!*************************************************************
! Regularisation function to be minimized for sqrt(s^2+x^2)
!        _      __________
!        \     /  2     2 |
! phi =  / \  / ss  +  x
!        -  \/           i
!        i
! (Carfantan: sqrtrf, with ss=0.01 or ss=0.001)
!*************************************************************
real function phi_sqrtrf(xx)
 real,dimension(idim)	:: xx

phi_sqrtrf = sum(sqrt(ss2 + xx**2))

return
end function phi_sqrtrf
!*************************************************************
! Gradient of potential function for convex sqrt(s^2+x^2) regularisation:
!        _      __________
!        \     /  2     2 |
! phi =  / \  / ss  +  x
!        -  \/           i
!        i
! (Carfantan: sqrtrf, with ss=0.01 or ss=0.001)
!
!*************************************************************
subroutine dphi_sqrtrf(xx,dx)
 real,dimension(idim)	:: xx,dx

 dx = xx/sqrt(ss2+xx**2)

 return
end subroutine dphi_sqrtrf
!*********************************************************
! Function E to be minimized for deconvolution 
!                                               
! This subroutine returns   E = |Y-Ax|^2 + alpha  phi
!                                              
!*********************************************************
real function func(xx)
! integer,parameter	:: idim=256
 real,dimension(idim)	:: xx
! The following statement makes the program crash
! real,dimension(:)	:: xx
 real			:: phi
 integer 		:: i

! Convolution = product in Fourier domain:
! w0 = convolution of xx with hh:
call dcv_conv_hh(xx,hh_re,hh_im,w0)
w0 = w0 - yy 

select case(option)
  case('tikho')
    phi = phi_tikho(xx) 
  case('maxen')
    phi = phi_mem(xx) 
  case('ggaus')
    phi = phi_ggauss(xx) 
  case('gmark')
    phi = phi_gmark(xx) 
  case('sqrtr')
    phi = phi_sqrtrf(xx) 
  case default
    print *,' func/Fatal error: invalid option'
    stop
end select

if(positivity)then
  do i = 1,nn
   if(xx(i) < 0.) phi = phi + 1000.*gscale
  end do
endif
func = sum(w0**2) + alpha * phi 

return
end function func
!*********************************************************
! Gradient of potential function for minimisation:
! grad_Er : value of the gradient of the criterium in X
!                                       (                 \---
! This subroutine returns  grad(E) =grad( |Y-Ax|^2+ alpha  >    phi(Xs-Xr)
!                                       (                 /___
!                                                      (s,r) in C
!
! d/dx (HX - Y)^T (HX - Y) = 2 H^T (HX - Y)
!*********************************************************
subroutine dfunc(xx,dx)
! integer,parameter	:: idim=256
 real,dimension(idim)	:: xx,dx
! The following statement makes the program crash
! real,dimension(:)	:: xx,dx

! Convolution = product in Fourier domain:
! w0 = conv(xx,hh):
call dcv_conv_hh(xx,hh_re,hh_im,w0)

dx = w0 - yy;

! w0 = conv(xx,hht):
call dcv_conv_hh(dx,hht_re,hht_im,w0)

select case(option)
  case('tikho')
    call dphi_tikho(xx,dx)
  case('maxen')
    call dphi_mem(xx,dx)
  case('ggaus')
    call dphi_ggauss(xx,dx)
  case('gmark')
    call dphi_gmark(xx,dx)
  case('sqrtr')
    call dphi_sqrtrf(xx,dx)
  case default
    print *,' dfunc/Fatal error: invalid option'
    stop
end select

if(positivity)then
  where(xx < 0.) dx = -1000. 
endif
 dx = 2*w0 + alpha*dx
 return
end subroutine dfunc
!*************************************************************
! Deconvolution by Wiener filter
!*************************************************************
subroutine dcv_wiener

! Compute FFT of signal: 
yy_re = yy; yy_im = 0.;
call fft_float(yy_re,yy_im,nx,ny,1);

! Compute transfer function:
hh_re = hh; hh_im = 0.;
call fft_float(hh_re,hh_im,nx,ny,1);

w0 = hh_re*hh_re + hh_im*hh_im;

! 1/(a+ib) = (a-ib)/(a2+b2)
! Wiener filter = ((a2+b2) / ((a2+b2) + alpha)) / (a+ib)
! Hence:        = (a-ib) / ((a2+b2) + alpha))
 w_re = hh_re / (w0 + alpha);
 w_im = - hh_im / (w0 + alpha);

! Deconvolution:
w0 = w_re*yy_re - w_im*yy_im
w1 = w_re*yy_im + w_im*yy_re
call fft_float(w0,w1,nx,ny,-1);
yyd = 0.; yyd(1:nn) = w0(1:nn)

return 
end subroutine dcv_wiener
!*************************************************************
! Deconvolution by spectral division 
!*************************************************************
subroutine dcv_spdiv

! Compute FFT of signal: 
yy_re = yy; yy_im = 0.;
call fft_float(yy_re,yy_im,nx,ny,1);

! Compute transfer function:
hh_re = hh; hh_im = 0.;
call fft_float(hh_re,hh_im,nx,ny,1);

w0 = hh_re*hh_re + hh_im*hh_im;

! 1/(a+ib) = (a-ib)/(a2+b2)
w_re = 0.; w_im = 0.;
where(w0 /= 0.)
 w_re = hh_re / w0;
 w_im = - hh_im / w0;
end where

! Deconvolution:
w0 = w_re*yy_re - w_im*yy_im
w1 = w_re*yy_im + w_im*yy_re
call fft_float(w0,w1,nx,ny,-1);
yyd = 0.; yyd(1:nn) = w0(1:nn)

return 
end subroutine dcv_spdiv
!*************************************************************
! Convolution H*X :  w0 = conv(xx,hh)
!
! IN:
! xx: function to be convolved by hh
! hh_re, hh_im: TF of hh
!
! OUT:
! ww = conv(hh,xx)  WARNING: often ww=w0 in the calling program...
!
!*************************************************************
subroutine dcv_conv_hh(xx,hh_re,hh_im,ww)
 real,dimension(idim)	:: xx,hh_re,hh_im,ww

! Convolution = product in Fourier domain:

! Compute FFT of xx: 
w_re = xx; w_im = 0.;
call fft_float(w_re,w_im,nx,ny,1);

! Product in Fourier domain:
ww = w_re*hh_re - w_im*hh_im
w1 = w_re*hh_im + w_im*hh_re

! Back to direct space:
call fft_float(ww,w1,nx,ny,-1);

return
end subroutine dcv_conv_hh
!*************************************************************
! Check the validity of the gradient
! and compare f(x+dx)-f(x)/dx with grad_f(x)
!
! x1,x2,dx: work space of dimension nn
!*************************************************************
subroutine dcv_check_grad(x1,x2,dx,nn)
 integer                :: nn,i
 real,parameter		:: eps=1.e-4,tolerance=1.e-4
 real			:: f_x2,f_x1,error
 real,dimension(idim)    :: x1,x2,dx

print *,' dcv_check_grad/Start checking gradient (Warning: may take some time...)'

! Generate random vector (between 0 and 1)
call random_seed; 
call random_number(x1);
f_x1 = func(x1) 

! Loop on all components:
!do i=1,nn
! Loop on a subsample (since too long otherwise):
do i=1000,1010
   x2 = x1;
   x2(i) = x1(i) + eps
   f_x2 = func(x1) 
   call dfunc(x1,dx)
   error = (f_x2 - f_x1)/eps - dx(i)*eps
   error = error / abs(f_x1 + 1.e-12)
   if(error > tolerance)then
      print *,'dcv_check_grad/Error!'
      write(6,12) i,error
12    format('i=',I4,' relative error =',E10.4)
   endif
end do

print *,' dcv_check_grad/End checking gradient'

return
end subroutine dcv_check_grad
!
!end module dcv_2D_mod
!***************************************************************
! module with frprmn
! Minimization of "func" with gradient "dfunc"
!
!***************************************************************
!module dcv_cgrad
! private
!! integer,parameter      	:: idim=256
!
! integer,private		:: ncom
! real,private,dimension(idim)	:: pcom,xicom
! public				:: frprmn
! private			:: norm_l2,linmin,df1dim,f1dim,mnbrak,dbrent
!
!contains
!*********************************************************
! Flechter-Reeves-Polak-Ribiere minimization
! Cf "Numerical recipees" Sect. 10.6 (in C p 423, in F77 p 416)
! IN:
! p: first guess
! n = dimension of p
! ftol = tolerance
! OUT:
! p: location of minimum
! iter: nber of iterations
! fret: value of minimum 
!
! Given a starting point p that is a vector of length n
! Flechter-Reeves-Polak-Ribiere minimization is performed on
! a function func, using its gradient as calculated by a routine dfunc.
! The convergence tolerance on the function value is input as ftol.
! Return quantities are p (the location of the minimum),
! iter (the number of iterations that were performed) 
! and fret (the minimum value of the function).
! The routine linmin is called to perform line minimizations. 
!********************************************************
      subroutine frprmn(p,n,ftol,iter,fret)
! itmax: maximum allowed number of iterations
! eps: small number to rectify the special case
! of converging to exactly zero function value
      integer itmax,j,iter,n
      real gam,gg,dgg,fret,ftol,eps,fp,norm_p,norm_p0
      real crit_fp,crit_p
      parameter (itmax=1000,eps=1.e-10)
      real p(n),g(idim),h(idim),xi(idim),p0(idim)
! Before
!      real func
!      external func
!      real norm_l2
!      external norm_l2
! Initializations:
      fp=func(p)
      call dfunc(p,xi)
      do 11 j=1,n
        g(j)=-xi(j)
        h(j)=g(j)
        xi(j)=h(j)
11    continue
      do 14 iter=1,itmax
        do j=1,n
          p0(j) = p(j)
        end do
        call linmin(p,xi,n,fret)
        do j=1,n
          p0(j) = p0(j) - p(j)
        end do
        norm_p0 = norm_l2(p0,n)
        norm_p = norm_l2(p,n)
        crit_p = norm_p0 / (norm_p+eps)
        crit_fp = 2.*abs(fret-fp)/(abs(fret)+abs(fp)+eps)
! For simulations every 10 iterations was OK, but not for real data...
!        if(mod(iter,10).eq.1)then
          write(6,18) iter,crit_p,crit_fp 
18	  format('Iter=',I3,' Rel. variation p=',e10.3,' fp=',e10.3)
!        endif
! Next statement is normal return:
! Exit criterium is 2 * |fret-fp| < ftol * (|fret| + |fp| + eps)
! i.e. relative variation of fp is small:
! 2 * |fret-fp| / (|fret| + |fp| + eps) < ftol
!  which is nearly: |fret-fp| / |fret| < ftol
!        if(2.*abs(fret-fp).le.ftol*(abs(fret)+abs(fp)+eps))return
! JLP 2001
! I also put a constraint on the variations of p
        if((crit_fp.le.ftol).and.(crit_p.le.ftol)) return
        fp=func(p)
        call dfunc(p,xi)
        gg=0.
        dgg=0.
        do 12 j=1,n
          gg=gg+g(j)**2
! Flechter-Reeves:
!         dgg=dgg+xi(j)**2
! Polak-Ribiere: (better)
          dgg=dgg+(xi(j)+g(j))*xi(j)
12      continue
! Unlikely: if gradient is exactly zero, then we are already done.
        if(gg.eq.0.)then
            write(6,*)'frprmn/return since gradient is null: gg=0'
            return
        endif
        gam=dgg/gg
        do 13 j=1,n
          g(j)=-xi(j)
          h(j)=g(j)+gam*h(j)
          xi(j)=h(j)
13      continue
14    continue
      write(6,*) 'frprmn/maximum iterations exceeded: iter=',iter
      return
      end subroutine frprmn
!************************************************************
! Norm L2
!************************************************************
real function norm_l2(p,n)
  integer 		:: n,j
  real,dimension(n)	:: p(n)
      norm_l2 = 0.
      do j=1,n
        norm_l2 = norm_l2 + p(j)*p(j)
      end do
      return
end function norm_l2
!************************************************************
! Search for minimum along a given direction
!************************************************************
subroutine linmin(p,xi,n,fret)
  integer		:: n,j
  real,parameter 	:: tol=1.e-4
  real,dimension(n)	:: p,xi
  real			:: fret,ax,xx,bx,xmin,fx,fa,fb
      ncom=n
      do j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
      end do 
      ax=0.
      xx=1.
      bx=2.
      call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
      fret=dbrent(ax,xx,bx,f1dim,df1dim,tol,xmin)
      do j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
      end do
      return
end subroutine linmin
!**************************************************************
real function df1dim(x)
  real,dimension(idim) 	:: xt,df
  integer		:: j
  real			:: x
      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
11    continue
      call dfunc(xt,df)
      df1dim=0.
      do 12 j=1,ncom
        df1dim=df1dim+df(j)*xicom(j)
12    continue
      return
end function df1dim
!**************************************************************
real function f1dim(x)
  real,dimension(idim) 	:: xt
  integer		:: j
  real			:: x
      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
11    continue
      f1dim=func(xt)
      return
      end function f1dim
!**************************************************************
subroutine mnbrak(ax,bx,cx,fa,fb,fc,func)
real,parameter 		:: gold=1.618034,glimit=100.,tiny=1.e-20
real			:: fa,fb,fc,func,fu,dum,ax,bx,cx
real			:: r,q,u,ulim
      fa=func(ax)
      fb=func(bx)
      if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+gold*(bx-ax)
      fc=func(cx)
1     if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),tiny),q-r))
        ulim=bx+glimit*(cx-bx)
        if((bx-u)*(u-cx).gt.0.)then
          fu=func(u)
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            go to 1
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            go to 1
          endif
          u=cx+gold*(cx-bx)
          fu=func(u)
        else if((cx-u)*(u-ulim).gt.0.)then
          fu=func(u)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+gold*(cx-bx)
            fb=fc
            fc=fu
            fu=func(u)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
          fu=func(u)
        else
          u=cx+gold*(cx-bx)
          fu=func(u)
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        go to 1
      endif
      return
end subroutine mnbrak
!**************************************************************
real function dbrent(ax,bx,cx,ff,dff,tol,xmin)
  integer,parameter	:: itmax=100
  real,parameter	:: zeps=1.0e-10
  logical		:: ok1,ok2
  integer		:: iter
  real			:: a,b,d,e,olde,u,v,w,x
  real			:: ax,bx,cx,ff,dff,tol,xmin,d1,d2,dv,dx,dw
  real			:: du,fu,fv,fx,fw,u1,u2
  real			:: xm,tol1,tol2
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
      fx=ff(x)
      fv=fx
      fw=fx
      dx=dff(x)
      dv=dx
      dw=dx
      do 11 iter=1,itmax
        xm=0.5*(a+b)
        tol1=tol*abs(x)+zeps
        tol2=2.*tol1
        if(abs(x-xm).le.(tol2-.5*(b-a))) goto 3
        if(abs(e).gt.tol1) then
          d1=2.*(b-a)
          d2=d1
          if(dw.ne.dx) d1=(w-x)*dx/(dx-dw)
          if(dv.ne.dx) d2=(v-x)*dx/(dx-dv)
          u1=x+d1
          u2=x+d2
          ok1=((a-u1)*(u1-b).gt.0.).and.(dx*d1.le.0.)
          ok2=((a-u2)*(u2-b).gt.0.).and.(dx*d2.le.0.)
          olde=e
          e=d
          if(.not.(ok1.or.ok2))then
            go to 1
          else if (ok1.and.ok2)then
            if(abs(d1).lt.abs(d2))then
              d=d1
            else
              d=d2
            endif
          else if (ok1)then
            d=d1
          else
            d=d2
          endif
          if(abs(d).gt.abs(0.5*olde))go to 1
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(dx.ge.0.) then
          e=a-x
        else
          e=b-x
        endif
        d=0.5*e
2       if(abs(d).ge.tol1) then
          u=x+d
          fu=ff(u)
        else
          u=x+sign(tol1,d)
          fu=ff(u)
          if(fu.gt.fx)go to 3
        endif
        du=dff(u)
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          dv=dw
          w=x
          fw=fx
          dw=dx
          x=u
          fx=fu
          dx=du
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            dv=dw
            w=u
            fw=fu
            dw=du
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
            dv=du
          endif
        endif
11    continue
      print *,'dbrent/Warning: exceeded maximum iterations.'
3     xmin=x
      dbrent=fx
      return
end function dbrent
!
!end module dcv_cgrad
end module dcv_2D_mod
