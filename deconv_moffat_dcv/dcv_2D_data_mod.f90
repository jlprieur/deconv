!****************************************************************** 
! Module used by dcv_deconv_2D.f90, dcv_cgrad_mod.f90 dcv_2D_mod.f90
!
! JLP
! Version 26/03/2001
!****************************************************************** 
module dcv_2D_data
public
! 16384=128*128
! integer,parameter		:: idim=128*128
 integer,parameter		:: idim=256*256
! integer,parameter		:: idim=512*512
 character(len=5)		:: option
!
! Used by dcv_deconv_2D et dcv_2D_mod:
!
 real,dimension(:),allocatable	:: hh,yy,yy0,yyd,yy_re,yy_im,w0,w1
 real,dimension(:),allocatable	:: hh_re,hh_im,hht_re,hht_im,w_re,w_im
 real				:: alpha
 integer			:: nx,ny,nn
 logical			:: positivity=.false.,simulation=.false.

end module dcv_2D_data
