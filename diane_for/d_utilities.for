C**************************************************************************
C Set of routines taken from SYLVIE.TLB
C Used by DIANE deconvolution routines
C Contains: FOURIER1, FOURIER, D_TO_S, S_TO_D, MOD_TF
C
C JLP
C Version 20-04-91
C--************************************************************************
	subroutine fourier1(br,bi,nx,ny,option)
C
C recoit l'image sous forme de deux tableaux R*8 de N valeurs :
C br est la partie reelle de l'image et bi la partie imaginaire (eventuelle)
C En sortie, la partie reelle et la partie imaginaire de la
C T.F. sont rangees de nouveau dans br et bi.
C
C OPTION :	+ 1	T.F. directe
C		- 1	T.F. inverse
C
	integer	nx,ny,npix(2),option, nsize, idim
 	real*8	br(nx,ny), bi(nx,ny)
	integer madrid(1),data
	common /vmr/madrid
C
	idim=nx
	npix(1)=nx
	npix(2)=ny
	nsize=8*nx*ny
C Allocation of memory:
	call jlp_getvm (data,nsize)
C Transfer:
	call move8tod(madrid(data),br,bi,nx,ny,idim)
C FFT
	call fourn1(madrid(data),npix,2,option)

C Transfer:
	call movedto8(madrid(data),br,bi,nx,ny,idim)

C Freeing memory:
	call jlp_freevm (data,nsize)
	return
        end
C********************************************************
C Subroutine to prepare the data for "fourn1"
	subroutine move8tod(data,br,bi,nx,ny,idim)
	real*4 data(*)
	real*8 br(idim,*),bi(idim,*)
	integer*4 i,j,k
	k=1
	do j=1,ny
	  do i=1,nx
          data(k)=sngl(br(i,j))
          data(k+1)=sngl(bi(i,j))
	  k=k+2
	  end do
	end do
        return
	end
C********************************************************
C Subroutine to sort out the data from "fourn1"
	subroutine movedto8(data,br,bi,nx,ny,idim)
	real*4 data(*)
	real*8 br(idim,*),bi(idim,*)
	integer*4 i,j,k
	k=1
	do j=1,ny
	  do i=1,nx
          br(i,j)=dble(data(k))
          bi(i,j)=dble(data(k+1))
	  k=k+2
	  end do
	end do
        return
	end
C*******************************************************
	subroutine fourier(br,bi,nx,ny,option)
C
C recoit l'image sous forme de deux tableaux R*8 de N valeurs :
C br est la partie reelle de l'image et bi la partie imaginaire (eventuelle)
C apres l'appel de C06FJF (NAG) la partie reelle et la partie imaginaire de la
C T.F. sont rangees de nouveau dans br et bi.
C
C OPTION :	+ 1	T.F. directe
C		- 1	T.F. inverse
C
	integer	dpwork,nx,ny,ntot,npix(2),option, wsize
 	real*8	br(nx,ny), bi(nx,ny)
	real*8	d_sn
C
	integer madrid(1)
	common /vmr/madrid
C
	ifail = 0
	ntot=nx*ny
        npix(1)=nx
        npix(2)=ny
C
C on calcule la racine de N en DP
C
	d_sn = sqrt (float(ntot))
C
C Pour le calcul de la T.F. inverse
C
	if (option .eq. -1) then
C Normalization for inverse FFT
	 do j = 1, ny	
	  do i = 1, nx
	   br(i,j) = br(i,j)/d_sn
	   bi(i,j) = bi(i,j)/d_sn
	  end do		
	 end do		
C Calling NAG routine:
CJLP92 	 call C06GCF (bi, ntot, ifail)
 	endif
C
C reservons de la place en memoire pour faciliter le calcul
C
	lwork = 3 * max(nx,ny)
        wsize = 8 * lwork
	call jlp_getvm (dpwork,wsize)
C
C calcul de la T.F. (NAG)
C (NAG divides by sqrt(ntot)...)
C
CJLP92	call C06FJF (2, npix, ntot, br, bi, madrid(dpwork), lwork, ifail)
C
C si on a calcule la T.F. inverse
C
	if (option .eq. -1) then
CJLP92	 call C06GCF (bi, ntot, ifail)
	else
C Normalization for direct FT
	 do j = 1, nx
	  do i = 1, ny     
	   br(i,j) = br(i,j)*d_sn
	   bi(i,j) = bi(i,j)*d_sn
	  end do
	 end do
	end if
C
	call jlp_freevm (dpwork,wsize)
	return
	end
C**************************************************************
	subroutine s_to_d (r4,r8,n)
C Passage de simple precision (R*4) en double precision (R*8)
C
	integer	n, i
	real*4		R4(n)
	real*8		R8(n)
C
	do i = 1,n
	 R8(i) = dble(R4(i))
	end do
	return
	end
C***********************************************************
	subroutine d_to_s (r8,r4,n)
C
C Passage de double precision (R*8) en simple precision (R*4)
C
	integer	n, i
	real*4	r4(n)
	real*8	r8(n)
C
	do i = 1,n
	 r4(i) = sngl(r8(i))
	end do
C
	return
	end
C***********************************************************
	subroutine mod_tf (n, preal, pimag)
C
C calcul du module
C
	integer	n, i
	real*8		preal(n), pimag(n)
C
	do i = 1,n
	 preal(i) = dsqrt((preal(i)*preal(i))+(pimag(i)*pimag(i)))
	end do
C
	return
	end
