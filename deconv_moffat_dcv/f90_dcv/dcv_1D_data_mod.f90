!****************************************************************** 
! Module used by dcv_deconv_1D.f90, dcv_cgrad_mod.f90 dcv_1D_mod.f90
!
! JLP
! Version 26/03/2001
!****************************************************************** 
module dcv_1D_data
public
 integer,parameter		:: idim=512
 character(len=5)		:: option
!
! Used by dcv_deconv_1D et dcv_1D_mod:
!
 real,dimension(:),allocatable	:: hh,yy,yy0,yyd,yy_re,yy_im,w0,w1
 real,dimension(:),allocatable	:: hh_re,hh_im,hht_re,hht_im,w_re,w_im
 real				:: alpha
 integer			:: nn,nn1,hh_nx

end module dcv_1D_data
