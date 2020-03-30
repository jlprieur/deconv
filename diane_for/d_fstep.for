*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*	PROGRAM D_FSTEP
*
* Central reconstruction with fixed steps with positivity constraint 
*
* INPUT:
*               First approx. of the object
*	*.V  	Support of the object (direct space)
*	*.G	Regularizing function
*
* OUTPUT:
*	*.PHIX	Final estimation of the object
*
* Log file: *.SUM
*
*
* Pointers :
*  ---------
* p_phit	phit		partie reelle			R4 -
* p_phitcr	phit chapeau	partie reelle			R4 -
* p_phitci	  "     "	partie imaginaire		R4 -
* p_v		support de l'Objet				R4 -
* p_g		Fonction de regularisation			R4 -
* p_phix	Objet reconstruit (module)			R4 -
* p_rn		Memory space for rn  				R4   vm
* p_dphitr	phit		partie reelle		R8	     vm
* p_dphiti	phit		partie imaginaire	R8	     vm
*
* From d_cgrad.for 
*
* JLP
* Version 02-05-91
*-------------------------------------------------------------------------------
	program d_fstep
	implicit none
	integer*4 npix(2),n,n4,n8,p_phit,p_phitcr, p_phitci
	integer*4 p_v,p_g,p_rn,p_phix,istatus
	integer*4 p_dphitr,p_dphiti,long,lu,i
	real*4 step
	character name*20,comments*80,filesum*40,buffer*90
 	character*40 fileout,filephit,filev,fileg
	integer*4 madrid(1)
	common /vmr/madrid
*
* Initialisation
	call jlp_begin
	call jlp_inquifmt
        write(6,112)
112     format('******* Program d_fstep  Version 02-11-93',
     1         ' Fixed step with positivity constraint')

*
* ouverture du fichier de la recapitulation
55	write(6,*)' Generic name (should end with a dot) ?'
	read(5,10) name
10	format(a)
	long=index(name,'.')-1
	if(long.le.0)goto 55
*
	filesum=name(:long)//'.SUM'
	lu=1
C "APPEND" is not available for IBM...
C	open(lu,file=filesum,status='OLD',access='APPEND')
	open(lu,file=filesum,status='OLD',access='SEQUENTIAL')
        do i=1,10000
        read(lu,10,end=989) buffer
        enddo
989	write (lu,980)
980	format(/,' MODULE FSTEP',/,' ----------------',/)
*
* Ouverture du fichier contenant la premiere approximation de l'Objet
*
	write(6,*)' First approximation of the object (test.PHIT?) :'
	read(5,10) filephit
	write(6,*)' Reading the first approx. of the object:',filephit
	call jlp_vm_readimag(p_phit,npix(1),npix(2),filephit,comments)
* Recentre image to have gravity center at (1,1):
        call recentre(madrid(p_phit),madrid(p_phit),
     1     npix(1),npix(2),npix(1))
*
	N = npix(1)*npix(2)
*
* reservons de la memoire et passons en D.P. pour les parties R et I
*
        N8 = N*8
	call jlp_getvm (p_dphitr,N8)
	call jlp_getvm (p_dphiti,N8)

* Memory space for FFT of first approximation:
        N4 = N*4
	call jlp_getvm (p_phitcr,N4)
	call jlp_getvm (p_phitci,N4)
*
* Ouverture du fichier contenant le support de l'objet
*
	FileV = name(:long)//'.V'
	write(6,*)' Reading the support of the object (direct space)',filev
	call jlp_vm_readimag(p_v,npix(1),npix(2),filev,comments)
* Recentre image to have gravity center at (1,1):
        call recentre(madrid(p_phit),madrid(p_phit),
     1     npix(1),npix(2),npix(1))
*
* Ouverture du fichier contenant la fonction de regularisation
*
	FileG = name(:long)//'.G'
	write(6,*)' Reading the regularization function (Fourier)',fileg
	call jlp_vm_readimag(p_g,npix(1),npix(2),fileg,comments)
* Recentre image to have zero frequency at (1,1):
        call recentre(madrid(p_g),madrid(p_g),
     1     npix(1),npix(2),npix(1))
*
	write(6,*) ' Step value for constant step gradients:'
	read(5,*) step 
*
* Memory working space: 
*
        n4 = n*4
        n8 = n*8
* Residuals:
	call jlp_getvm (p_rn,n4)
* FFT of first approximation:
	call jlp_getvm (p_phitcr,N4)
	call jlp_getvm (p_phitci,N4)
* Final approximation:
	call jlp_getvm(p_phix,N4)
 
*
C	   call restit_cgrad (N,npix,madrid(p_phitcr),madrid(p_phitci),
C     1	madrid(p_phit),madrid(p_dphitr),madrid(p_dphiti),
C     1	madrid(p_v),madrid(p_g),madrid(p_rn),madrid(p_zn),
C     1	madrid(p_ddnr),madrid(p_ddni),madrid(p_ddntemp),
C     1	madrid(p_phix),p_dphitr,p_dphiti,
C     1	p_ddnr,p_ddni,eps,lu)
* Fixed step gradient with positivity constraint:
	   call restit_fstep (N,npix,madrid(p_phitcr),madrid(p_phitci),
     1	madrid(p_phit),madrid(p_dphitr),madrid(p_dphiti),
     1	madrid(p_v),madrid(p_g),madrid(p_rn),
     1	madrid(p_phix),lu,step)
*
* et on ferme les fichiers
*
	fileout = name(:long)//'.PHIX'
	comments = 'Resultat final'
C* Recentre image to have gravity center at (1,1):
C        call recentre(madrid(p_phix),madrid(p_phix),
C     1     npix(1),npix(2),npix(1))
	call jlp_writeimag(madrid(p_phix),npix(1),npix(2),npix(1),
     1	fileout,comments)
*
* Sortie normale
	close(lu)
	call jlp_end
	stop
	end
*------------------------------------------------------------------------------
*
*-----------------------------------------------------------------------
* Like restit_cgrad but with fixed step and positivity 
* constraint
* JLP May 1992
*
*-----------------------------------------------------------------------
	subroutine restit_fstep(n,npix,phitcr,phitci,phit,
     1	dphitr,dphiti,v,g,rn,phix,lu,step)
	implicit none
       	integer*4  n,npix(2),lu,imax
	integer*4 k,i,m,mult(35),icn(101),isig,num,nomb,nrac
        real*4 epsilon
	real*4 v(n),g(n),rn(n),phix(n),phitcr(n), phitci(n),phit(n)
 	real*4 rac(35),xnorfin,xx,x,xnun,pas,xnfinc
        real*4 deltan,ron,tete,tetap,xmum,step
	real*8 dphitr(n), dphiti(n)
	real*8 dg2,xnorrn,psdz,omegan
	real*8 error,major,y,xnorrnp	

* Majorant de 1/mu:
	major=10.	
*
*
*
* Computing FFT of the first approximation of the object:
*
* Transfer to work arrays:
	do k = 1, n
	 phitcr(k)=phit(k)
	 phitci(k)=0.
	end do
* Compute phi_n chap
        ioption=1
	call fft_2d(phitcr,phitci,npix(1),npix(2),npix(1),ioption)
************************************************************************
*
*
* N = npix ** 2       nombre total de pixels
* (phitcr,phitci)     tf de phit - parties reelle et imaginaire
* (phit)	     module de phit
* dphitr et dphiti  tableaux associes en dp
* V		   support de l'objet
* G		   fonction de regularisation
*
************************************************************************
*
* initialisation
* of first guess  phi_zero = v * phit
	do k = 1, n
	 phix(k)=v(k)*phit(k)	
	end do
*******************************************************************************
*
*			I T E R A T I O N
*
*******************************************************************************
*
	write(6,105)
	write(lu,105)
105	format('  N         Error',/)
*
	imax = 30 
        epsilon = 0.1
        do 600 iter=1,imax
*
* Compute residuals: r_n = || aa*(psi - aa phi)||
*  i.e., r_n = modulus of ( v U* g2 (phitc - phi_n chap))
*

* Transfer to work arrays:
	xnorrn = 0.
	do k = 1, n
	 dphitr(k)=dble(phix(k))
	 dphiti(k)=0.
	 xnorrn = xnorrn + dphitr(k) * dphitr(k)
	end do
C        print *,' norm of FFT of phit troncated:',xnorrn

* Compute phi_n chap
        ioption=1
	call fourier1(dphitr,dphiti,npix(1),npix(2),ioption)
*
	xnorrn = 0.
        work = 0.
	do k = 1, n
	 dg2 = dble(g(k)) * dble(g(k))
	 dphitr(k) = (dble(phitcr(k)) - dphitr(k)) * dg2
	 dphiti(k) = (dble(phitci(k)) - dphiti(k)) * dg2
	 xnorrn = xnorrn + dphitr(k) * dphitr(k)
         work = work + dg2
	end do
C        print *,' norm of FFT of initial phit - FFT of phit troncated:',xnorrn
C        print *,' norm of dg2:',work
* Inverse FFT:
        ioption=-1
	call fourier1(dphitr,dphiti,npix(1),npix(2),ioption)

* Compute FFT modulus, and store it in real part array:
* I assume zero phase (real image), so I take the modulus of the inverse FFT:
	call mod_tf(n,dphitr,dphiti)

* Compute residuals:
	xnorrn = 0.
	do k = 1, n
	 rn(k) = sngl(dphitr(k)) * v(k)
	 xnorrn = xnorrn + rn(k) * rn(k)
	end do
*
*
* Compute phi_n+1 and ||phi_n+1||**2
* Compute r_n+1 and ||r_n+1||**2
*
	xnorfin = 0.
*
	do k = 1, N
C New approximation:
	 phix(k) = phix(k) + step * rn(k)
C Positivity constraint:
         if(phix(k).lt.0.)phix(k)=0.
	 x = phix(k) * phix(k)
	 xnorfin = xnorfin + x
	end do
*
* test d'arret
*
	error = xnorrn/xnorfin
	error = major * sqrt(error)
*
	write (6,102) iter, error
	write (lu,102) iter, error
102	format(1x, i3, 4x, g12.6)
*
* Test to exit from the loop: 
	if (error .le. epsilon) goto 500

* Go to beginning of loop:
600     continue
*
* ********************  FIN D'ITERATION   ********************
*
        write(6,28)
        write(lu,28)
28      format('jlp_fstep/ Warning: error larger than expected!')

500	return
	end
C-----------------------------------------------------------------------
	subroutine dnormel2(n,reel,imagi,norme)
*
* Calcul de la norme L2 du complexe (REEL, IMAG) dans l'espace de Fourier
*
* NE DOIT PAS ETRE APPELE EN UTILISANT LE MECANISME	%VAL
* ----------------------------------------------------------
*
	implicit none
	integer*4 n,i
	real*8	reel(*), imagi(*)
	real*4	norme
*
	norme = 0.
	do i = 1, n
	 norme = norme + (reel(i)*reel(i)+imagi(i)*imagi(i))
	end do
*
	norme = sqrt(norme/float(n))
*
	return
	end
C-----------------------------------------------------------------------------
C	include 'd_utilities.for'
