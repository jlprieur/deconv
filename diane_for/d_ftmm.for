*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*	PROGRAM D_FTMM
*
* Normalisation L1 de la reponse impulsionnelle
* Transformee de Fourier de la reponse impulsionnelle normalisee
*
* New version: input is square modulus and zero frequency is in the center
*
* Input: square modulus of the Transfer Function 
* Output: *.FTB, modulus of the Bounded Transfer Function
*
* From d_ftm S. ROQUES - J. VIGNEAU	30 NOVEMBRE 1989
*
* JLP (Modif Midas)
* Version 12-11-96
*-------------------------------------------------------------------------------
	program d_ftmm
	integer*4 nx,ny,pntrin,pntrout
	integer*4 isize,i,long,lu
	real*4 thresh	
	character comments*80,buffer*90
	character*40 filein,fileout,filesum
	integer*4 madrid(1)
	common /vmr/madrid
*
	call jlp_begin
*
* Initialisation
*
	call jlp_inquifmt
	write(6,102)
102     format('******* Program d_ftmm  Version 12-11-96')

	write(6,12)
12	format(' Normalization of the Bounded Transfer Function',/,
     1         ' Input: modsq_* , square modulus of the Transfer Function',/,
     1         '(zero frequency is at the central pixel)',/,
     1         ' Output: *.FTB, modulus of the Bounded Transfer Function')
*
* Lower threshold for PSF support in Fourier domain
	write(6,*) ' Lower threshold of PSF in Fourier domain (0.01?)' 
	read(5,*) thresh 
*
* on ouvre le fichier recapitulation
*
	lu=1
C	FileSum= (name(:long)//'.SUM')
	FileSum='ftmm.log' 
C "APPEND" is not available for IBM...
CC       open(lu,file=filesum,status='OLD',access='APPEND')
C        write(6,28) FileSum
C28      format(' Opening (old) log file: ',a,/,
C     1         ' In case of error, please create such a file')
C        open(lu,file=filesum,status='OLD',access='SEQUENTIAL')
        open(lu,file=filesum,status='unknown',access='sequential')
        do i=1,10000
        read(lu,10,end=989) buffer
10      format(a)
        enddo
989	write (lu,980)
980	format(/,' Program d_ftmm (version 12-11-96)',/,' ----------')
*
* Ouverture du fichier contenant la reponse impulsionnelle
*
	write(6,*) ' Input square modulus of the PSF (Fourier space) :'
	read(5,10) filein
	call jlp_vm_readimag(pntrin,nx,ny,filein,comments)
	write (lu, '(a)') (' Square modulus of the PSF : '//filein)
	write (lu, '(a)') (' '//comments)
*
* on normalise la reponse impulsionnelle
*
	isize = nx * ny
	i=isize*4
	call jlp_getvm(pntrout,i)
	call sqroot(nx,ny,nx,madrid(pntrin),madrid(pntrout))
*
* Apply a threshold constraint on the transfer function 
*
	call borne(isize,madrid(pntrout),thresh)
*
* and output of the Bounded Transfer Function
*
C	fileout = (name(:long)//'.FTB')
	write(6,*)' Writing the Bounded Transfer Function (modulus):'
	read(5,10) fileout
	write(comments,11) filein(1:20),thresh
11      format('Bounded Tr.F. ',A20,'/thresh:',G10.3)
	call jlp_writeimag(madrid(pntrout),nx,ny,nx,fileout,comments)
*
* Sortie normale
*
	close(lu)
	call jlp_end
	stop
	end
*
*-------------------------------------------------------------------------------
* Apply a threshold constraint on the transfer function 
*
	subroutine borne(n,a,thresh)
	implicit none
	integer*4 n
	real*4 	a(n), thresh 
	integer*4 i
*
	do i = 1,n
	 if(a(i).le.thresh) a(i)=0.
 	end do
*
	return
	end
*
C-----------------------------------------------------------------------------
C Computes square root of square modulus
C-----------------------------------------------------------------------------
	subroutine sqroot(nx,ny,idim,input,output)
	real input(idim,*),output(idim,*)
	real w1
	integer nx,ny,icenter,jcenter

C Assumes zero frequency is in the center (nx/2+1,ny/2+1) 
        icenter=nx/2+1
        jcenter=ny/2+1
	w1=input(icenter,jcenter)
	print *,' Initial value of zero frequency:',w1
	do j=1,ny
	  do i=1,nx
	   w2=input(i,j)/w1
	   if(w2.gt.0.)then
	     output(i,j)=sqrt(w2)
	   else
	     output(i,j)=0.
	   endif
	  end do
	end do
	return
	end
C-----------------------------------------------------------------------------
C	include 'd_utilities.for'
