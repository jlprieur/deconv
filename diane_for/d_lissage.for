*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*	PROGRAM D_LISSAGE
*
* Compute smoothing function and raw approximation of the object.
* From lissage / S. ROQUES - J. VIGNEAU		30 NOVEMBRE 1989
*
* Input/output FFT: zero frequency at the center.
*
* INPUT:
*	*.FTB   Transf. function (modulus)
*	*.SNR	SNR of the input image (Fourier domain)
*	*.RFT	FT of the input image (real part)
*	*.IFT	FT of the input image (imag. part)
*
* OUTPUT:
*	*.HR	Synthetic aperture
*	*.FLI   Smoothing filter
*	*.PHIT  First approx. of the object
*	*.PCR   First approx. of the object: real part of FT
*	*.PCI   First approx. of the object: imag. part of FT
* 
* JLP
* Version 02-10-96
*-------------------------------------------------------------------------------
	PROGRAM D_LISSAGE
	implicit none
*
	integer*4 npix(2), npix1(2), npix2(2),
     1	pntrs,pntrhr,pelipr,pftm,psnr,prft, pift,
     1	phi_tcr,phi_tci,phi_tr,dphi_tr,dphi_ti,
     1	dpstock,dplisse,dptempi
	integer*4 n,n4,n8,i,lu,long,istatus
	real*4	fradius,sector_angle_deg,alphat,Nor_phi_T
	character*40 FileIn,FileOut,FileOut2,FileOut3,FileOut4,FileHr
	character name*20,rep*1,comments*80,FileSum*40,buffer*90
	integer*4 madrid(1)
	common /vmr/madrid
*
* Initialisation
*
	call jlp_begin
	call jlp_inquifmt
        write(6,102)
102     format('******* Program d_lissage  Version 02-10-96')
*
* saisie de l'angle (cas de la radioastronomie)
*
C	write(6,12)
C12	format(' Angle of the sector (for radio-astronomy,',
C     1	' otherwise 180.) ?')
C	read(5,*) sector_angle_deg
        sector_angle_deg=180.
*
* saisie de la valeur seuil pour snr (alphaT)
*
	write(6,*) ' SNR threshold to reject the spectral points (alphaT) ?'
	read(5,*) alphat
* Storing this value in the symbol D_ALPHAT
	write(buffer,*) alphat
	call jlp_set_symbol('D_ALPHAT',buffer,istatus)
*
* Generic name:
55	write(6,*) ' Generic name ? (end with a dot)'
	read(5,10) name
10	format(a)
	long=index(name,'.')-1
	if(long.le.0) goto 55
* Radius of the support Hr in Fourier domain:
	write(6,*) ' Radius you want to reach in Fourier domain (radius of Hr)?'
	read(5,*) fradius 
*
* on ouvre le fichier de la recapitulation
*
	FileSum=name(:long)//'.SUM'
	lu=1
C "APPEND" is not available for IBM...
C       open(lu,file=filesum,status='OLD',access='APPEND')
        open(lu,file=filesum,status='OLD',access='SEQUENTIAL',err=879)
881     do i=1,10000
        read(lu,10,end=989) buffer
        enddo
989	write (lu,980)
980	format(/,' Program d_lissage',/,' --------------')
*
	write (lu, 1002) fradius 
1002	format (' Radius of the support Hr in Fourier space : ', f5.2)
	write (lu, 1003) alphat
1003	format (' Threshold for SNR in Fourier space (alphaT) : ', f5.2)
C	write (lu, 1004) sector_angle_deg
C1004	format (' Angle : (deg)', f6.2)
*
* ouverture du fichier de la fonction de Transfert en F_FTB
	write(6,*) ' Reading the limited Transfer function (*.FTB):'
	read(5,10) FileIn
C        FileIn=name(:long)//'.FTB'
	call jlp_vm_readimag(pftm,npix(1),npix(2),FileIn,comments)

C Recentering the image: to get zero frequency at (1,1):
        call recentre(madrid(pftm),madrid(pftm),npix(1),npix(2),npix(1))

* Precaution avant de commencer
* Power of two only!
        iok_flag = 0
        jok_flag = 0
        do i=1,12
           if(npix(1).eq.2**i)iok_flag =1
           if(npix(2).eq.2**i)jok_flag =1
        enddo
        iok_flag = iok_flag * jok_flag
	if (iok_flag.ne.1) then
	 write(6,*) 'Fatal error/Only powers of two are allowed for nx and ny!'
	 stop 
	end if
C	if (npix(1).ne.npix(2)) then
C	 write(6,*) 'On ne peut travailler que sur un carre'
C	 stop 
C	end if
*
*
***************************************************************
* CALL LISSE:
*
* Reservons de la place
	N=npix(1)*npix(2)
	n4=n*4
	n8=n*8
*
* Getting memory space for smoothing filter
	call jlp_getvm(pntrs,n4)
* Getting memory space for synthetic aperture Hr:
	call jlp_getvm(pntrhr,n4)
	call jlp_getvm(pelipr,n4)
*
	call jlp_getvm (dpstock,n8)
	call jlp_getvm (dplisse,n8)
	call jlp_getvm (dptempi,n8)
*
* Compute dplisse, smoothing function in Fourier domain: 
*
	call lisse(fradius,npix(1),npix(2),sector_angle_deg,madrid(pntrhr),
     1   madrid(pelipr),madrid(dpstock),madrid(dplisse),
     1    madrid(dptempi),dplisse, dptempi,lu)
*
* Conversion of "dplisse" from double to simple precision 
* (Output is "pntrs")
*
	call D_to_S (madrid(dplisse), madrid(pntrs), N)
*
	call jlp_freevm (pelipr,n4)
	call jlp_freevm (dpstock,n8)
 	call jlp_freevm (dptempi,n8)
 	call jlp_freevm (dplisse,n8)
*
* Synthetic aperture:
	FileHr = name(:long)//'.HR'
	write(6,*) ' Writing the synthetic aperture: ',FileHr
	comments = 'Synthetic aperture Hr (Fourier) for '//name(:long)
C Recentering the image: to get zero frequency in the center:
        call recentre(madrid(pntrhr),madrid(pntrhr),npix(1),npix(2),npix(1))
	call jlp_writeimag(madrid(pntrhr),npix(1),npix(2),
     1	npix(1),filehr,comments)
 	call jlp_freevm (pntrhr,n4)

*
*********************************************************************
*
* calcul de phit
*
*	POINTEURS :
*	---------
*	pftm	fonction de transfert	<name>.FTB
*	psnr
*	prft
*	pift
*	pntrs	fonction de lissage en R*4
*
* Fichiers de sortie
*	- phi_T_chapeau		partie reelle
*	-	"		partie imaginaire
**************************************************************************
*
* Ouverture du fichier contenant le rapport signal/bruit
*
	write(6,*) ' Reading the SNR of the image (Fourier domain) *.SNR'
 	FileIn = name(:long)//'.SNR'
	call jlp_vm_readimag(psnr,npix1(1),npix1(2),FileIn,comments)
        write(6,*) " nx, ny ",npix1(1),npix1(2)
C Recentering the image: to get zero frequency at (1,1):
        call recentre(madrid(psnr),madrid(psnr),npix1(1),
     1   npix1(2),npix1(1))
*
* ouverture du fichier contenant la partie reelle de la T.F. de l'image
*
	write(6,*) ' Reading the real part of F.T. of the image *.RFT'
 	FileIn = name(:long)//'.RFT'
	call jlp_vm_readimag(prft,npix2(1),npix2(2),FileIn,comments)
C Recentering the image: to get zero frequency at (1,1):
        call recentre(madrid(prft),madrid(prft),npix2(1),
     1   npix2(2),npix2(1))
*
* ouverture du fichier contenant la partie imaginaire...
*
	write(6,*) ' Reading the real part of F.T. of the image *.IFT'
 	FileIn = name(:long)//'.IFT'
	call jlp_vm_readimag(pift,npix2(1),npix2(2),FileIn,comments)
C Recentering the image: to get zero frequency at (1,1):
        call recentre(madrid(pift),madrid(pift),npix2(1),
     1    npix2(2),npix2(1))

*
* on verifie que les images sont compatibles
*
	do i = 1, 2
 	 if (npix(i) .ne. npix1(i).or.
	1    npix(i) .ne. npix2(i)) then
	     write(6,13)
13	     format(' d_lissage/Fatal error: ',/,
     1	' Input images have not the same size!')
	     stop 
	 endif
	end do
*
* Premiere approximation de l''objet - T.F. partie reelle (FileOut2):
	call jlp_getvm(phi_tcr,n4)
*
* Premiere approximation de l''objet - T.F. partie imaginaire (FileOut3):
	call jlp_getvm(phi_tci,n4)
*
* Premiere approximation de l'objet FileOut4:
	call jlp_getvm(phi_tr,n4)
*
	call jlp_getvm (dphi_tr,n8)
	call jlp_getvm (dphi_ti,n8)
* phi_t(n,alphat,s,h,snr,phicr,phici,phitcr,phitci,dphitr,
*        dphiti,nor_phi_t,phitr)
*
 	call phi_t(n,alphat,madrid(pntrs),madrid(pftm),madrid(psnr),
     1	    madrid(prft),madrid(pift),madrid(phi_tcr),madrid(phi_tci),
     1	    madrid(dphi_tr),madrid(dphi_ti),nor_phi_t,madrid(phi_tr))
*
	write (6,2005) nor_phi_t
	write (lu,2005) nor_phi_t
2005	format(' NORM of smoothing function = ',F12.4)

*
* Smoothing filter:
	fileout = name(:long)//'.FLI'
	write(6,2003) fileout
2003	format(/,' Writing the smoothing filter (Fourier domain):',A)
	write (comments,1000) fradius 
1000	format (' Smoothing filter (Fourier domain) - Hr radius = ',F6.2)
C Recentering the image: to get zero frequency in the center:
        call recentre(madrid(pntrs),madrid(pntrs),npix(1),npix(2),npix(1))
	call jlp_writeimag(madrid(pntrs),npix(1),npix(2),
     1	npix(1),FileOut,comments)
 
* Premiere approximation de l''objet - T.F. partie reelle (FileOut2):
	FileOut2 = name(:long)//'.PCR'
	write(6,2000) FileOut2
2000	format(/,' Writing the first approx. of the F.T. of the object',
     1	/,' Real part :',A)
	comments = 'First approximation of the FT of the object (real)'
C Recentering the image: to get zero frequency in the center:
        call recentre(madrid(phi_tcr),madrid(phi_tcr),
     1    npix(1),npix(2),npix(1))
	call jlp_writeimag(madrid(phi_tcr),npix(1),npix(2),
     1	  npix(1),FileOut2,comments)
*
* Premiere approximation de l''objet - T.F. partie imaginaire (FileOut3):
	FileOut3 = name(:long)//'.PCI'
	write(6,2001) FileOut3
2001	format(/,' Writing the first approx. of the F.T. of the object',
     1	/,' IMAGinary part :',A)
	comments = 'First approximation of the FT of the object (imag)'
C Recentering the image: to get zero frequency in the center:
        call recentre(madrid(phi_tci),madrid(phi_tci),
     1    npix(1),npix(2),npix(1))
	call jlp_writeimag(madrid(phi_tci),npix(1),npix(2),
     1	  npix(1),FileOut3,comments)
 
* First approximation of the object:
	FileOut4 = name(:long)//'.PHIT'
	write(6,2002) FileOut4
2002	format(/,' Writing the first approx. of the object :',A)
	comments = ' First approximation of the objet'
C Recentering the image: to get zero frequency in the center:
        call recentre(madrid(phi_tr),madrid(phi_tr),npix(1),npix(2),npix(1))
	call jlp_writeimag(madrid(phi_tr),npix(1),npix(2),
     1	  npix(1),FileOut4,comments)
*
* Sortie normale
	close(lu)
	call jlp_end
	stop
* Error (logfile has not been created yet):
879     write(6,878) FileSum 
878     format(' >>>>>>> D_LISSAGE/Error opening (old) logfile ',A)
        open(lu,file=filesum,status='NEW',access='SEQUENTIAL')
        goto 881

	end
*-------------------------------------------------------------------------------
* lisse computes the smoothing function "tempr"
* with successive FFT and truncation both in direct and Fourier space
*
* INPUT:
* fradius: limit in resolution (Hr radius in Fourier space) 
*
*-------------------------------------------------------------------------------
 	subroutine lisse(fradius,nx,ny,sector_angle_deg,sectr,elipr,
     1	stock,tempr,tempi,dptempr,dptempi,lu)
	implicit none
*
	integer*4 nx,ny,dptempr, dptempi,lu
	real*4 sectr (nx,*),elipr (nx,*),fradius,Gain,sector_angle_deg
	real*8 stock (nx,*),tempr (nx,*),tempi (nx,*)
	integer*4 N2,ixcenter,iycenter,n_sur_2,k1,l1
        integer*4 irp,r2p,is, ic,l, k,c2, d2,test_Int,itere, it_max
 	real*4 sector_angle_rad,test_real,x,A,eta0,pi,rr,Rp
        real*4 moment_Y, moment_X,rap,gr_axe, p_axe
	real*8 error,epsilon,norme,X2,Y2
	integer*4 madrid(1)
	common /vmr/madrid

C eta0 = 1.55 for Khi**2=0.95
	eta0=1.55
	pi=3.14159265358979323846
	it_max=50
	epsilon=1.E-5
*
	n2 = nx * ny
	ixcenter = nx/2 + 1
	iycenter = ny/2 + 1
C limite de resolution (fradius - espace de Fourier) et rayon du disque
C (rr - espace reel) si la configuration etait disque-disque
C JLP92: lr = eta0*nx/(rr*pi) with eta0 = 1.55
C Comment: eta0 = rr * (fradius / nx) * pi
	 rr = eta0*nx/(fradius*pi)
	 write (6, 1002) rr
	 write (lu, 1002) rr
1002	format (' I set the limit of resolution so that eta0=1.55 (disks)',/,
     1         ' We find here that the radius of the smoothing function',
     1         ' (direct space) is: ',F6.3,' pixels')
C JLP93:
C If fradius in Fourier space 
C (i.e. fradius/nx spatial freq. thus radius = nx/fradius in direct space) 

*
* calcul du rayon du secteur donnant la meme surface que le disque
*
	rp = fradius * sqrt (180./sector_angle_deg)
	irp = int(rp) + 1
	r2p = irp * irp
	is = 0	
* Conversion to radians:
	sector_angle_rad = sector_angle_deg * pi/180.
*
* Computation of the Fourier "mask",
* i.e. the characteristic function (with sectors), SECTR(K,L)
*
	if (sector_angle_deg.ne.180.) test_real = tan(sector_angle_rad/2.)
*
   	do  j = 1, ny
	 c2 = (j-iycenter) * (j-iycenter)
	 do i = 1, nx
	  sectr (i,j) = 0.
	  test_int = c2 + (i-ixcenter)**2
	  if (test_int .le. r2p) then
	   if (sector_angle_deg .eq. 180.) then
* cas du disque
	    sectr (i,j) = 1.
	    is = is + 1
	   else			
* cas du secteur
	    if (i .eq. ixcenter) then
	     sectr (i,j) = 1.
	     is = is + 1
	    else
	     if (j .eq. iycenter) then
	      A = 0.
	     else
	      A = abs(float(j-iycenter)/float(i-ixcenter))
	     endif
	     if (a .le. test_real) then
	      sectr (i,j) = 1.
	      is = is + 1
	     endif
	    endif
	   endif
	  endif
	 end do
 	end do
*
* on recentre pour travailler dans l'espace des frequences
*
	n_sur_2 = nx/2
	do j = 1, n_sur_2
	 j1 = j + n_sur_2
	 do i = 1, n_sur_2
	  i1 = i + n_sur_2
	  x = sectr (i1, j1)
	  sectr (i1, j1) = sectr (i, j)
	  sectr (i, j) = x
	  x = sectr (i1, j)
	  sectr (i1, j) = sectr (i, j1)
	  sectr (i, j1) = x
	 end do
	end do
*
* calcul de l'ellipsoide de resolution dans l'espace reel
* et fonction caracteristique, ELIPR(K,L)
*
	moment_Y = sector_angle_rad + sin(sector_angle_rad)
	moment_X = sector_angle_rad - sin(sector_angle_rad)
	rap = moment_Y/moment_X
	p_axe = rr
	gr_axe = p_axe * rap
	ic = 0
*
	p_axe = p_axe*p_axe
	gr_axe = gr_axe*gr_axe
	do j = 1, ny
	 c2 = (j - iycenter) * (j - iycenter)
	 do i = 1, nx
 	  elipr(i,j) = 0.
	  d2 = (i-ixcenter) * (i-ixcenter)
	  test_real = float(c2)/p_axe + float(d2)/gr_axe
	  if (test_real .le. 1.) then
	   elipr(i,j) = 1.
	   ic = ic + 1
	  endif
	  tempr(i,j) = dble (elipr(i,j))
	 end do
	end do
*
* reajustement de eta0
*
	eta0 = (sqrt(float(ic*is)))/nx
*
	write (6,999) eta0
	write (lu,999) eta0
999	format (' Actual value of interpolation coefficient eta0: ',F6.3,
     1          ' (Do not worry as long as eta0 is within [1.4 , 2.2])')
	if (p_axe .ne. gr_axe) then
	 write (6,1000) sqrt(p_axe)
	 write (lu,1000) sqrt(p_axe)
1000	format (' Petit axe de la fonction de lissage : ',F6.3)
	 write (6,1001) sqrt(gr_axe)
	 write (lu,1001) sqrt(gr_axe)
1001	format (' Grand axe de la fonction de lissage : ',F6.3)
	else
	 write (lu, '(a)') ' '
	end if
	write (6, 1005) fradius 
	write (lu, 1005) fradius 
1005	format (' Radius of Hr (Fourier space)  : ',F8.3)
	write (lu, '(a)') ' '
*
* Calcul de la fonction de lissage par la methode de la puissance
*
	itere = 0
	error = epsilon + 1.
*
C JLP91	do while (itere .lt. it_max .and. error .gt. epsilon)
*
88       continue
*
	 do j = 1, ny
	  do i = 1, nx
	   stock(i,j) = tempr(i,j)
	   tempi(i,j) = 0.
	  end do
	 end do
*
*... T.F. and truncation in Fourier space
*
* JLP96 I remove
*	 call fourier1(madrid(dptempr),madrid(dptempi),nx,ny,1)
	 call fourier1(tempr,tempi,nx,ny,1)
*
	 do 303 j = 1, ny
	  do 302 i = 1, nx
 	   if (sectr(i,j) .eq. 0.) then
	    tempr(i,j) = 0.
	    tempi(i,j) = 0.
	   end if
302	  continue 
303	 continue 
*
*... Inverse FFT and truncation in direct space:
*
* JLP96 I remove
*	 call fourier1(madrid(dptempr),madrid(dptempi),nx,ny,-1)
	 call fourier1(tempr,tempi,nx,ny,-1)
*
	 do j = 1, ny
	  do i = 1, nx
	   if (elipr(i,j) .eq. 0.) then
	    tempr(i,j) = 0.
	    tempi(i,j) = 0.
	   end if
	  end do
	 end do
*
*... normalisation et test de convergence
*
	 norme = 0.
	 error = 0.
*
	 do j = 1, ny
	  do i = 1, nx
	   X2 = tempr(i,j) * tempr(i,j)
	   Y2 = tempi(i,j) * tempi(i,j)
	   norme = norme + X2 + Y2
	   tempr(i,j) = DSQRT (X2 + Y2)
	  end do
	 end do
*
	 norme = DSQRT (norme)
*
	 do l = 1, ny
	  do k = 1, nx
	   tempr (k,l) = tempr (k,l)/norme
	   X2 = (tempr(k,l) - stock(k,l)) * (tempr(k,l) - stock(k,l))
	   error = error + X2
	  end do
	 end do
*
	 error = DSQRT (error)
 	 itere = itere + 1
*
C JLP91
C JLP91	do while (itere .lt. it_max .and. error .gt. epsilon)
C	end do
         if (itere .lt. it_max .and. error .gt. epsilon) goto 88
********************************************************************
* End of loop:
*
* Normalization (L1) in direct space:
	do j = 1, ny
	 do i = 1, nx
	  tempi (i,j) = 0.
	 end do
	end do
	call dnormel1(nx,ny,madrid(dptempr),madrid(dptempr))
*
* Final smoothing function is taken as the modulus of the FFT of that
* normalized function:
	call fourier1(madrid(dptempr),madrid(dptempi),nx,ny,1)
	call mod_tf(n2, madrid(dptempr), madrid(dptempi))
*
	write (6, 1003) itere,it_max
	write (lu, 1003) itere,it_max
1003	format(' Number of iterations: (computation of smoothing function)',
     1          i3,' (it_max=',i3,')')
	write (6, 1004) error, norme
	write (lu, 1004) error, norme
1004	format (' Erreur : ', e8.3, '              Khi  : ', f5.3)
*
	return
	end
*
*-------------------------------------------------------------------------------
*
 	subroutine phi_t(n,alphat,s,h,snr,phicr,phici,
     1	phitcr,phitci,dphitr,dphiti,nor_phi_t,phitr)
*
* Premiere approximation de l'objet.
*
* N	: Nombre total de pixels				Integer
* alphat: Seuil pour le rapport Signal/Bruit			real
* S	: Fonction de lissage dans l'espace des frequences	Tableau reel
* H	: Fonction de transfert		"	"		    "
* snr	: Rapport Signal/Bruit		"	"		    "
* phicr	: TF de l'image experimentale (partie reelle)	 	    "
* phici :	"	"	"     (partie imaginaire)	    "
* phitcr: TF de la premiere approximation (partie reelle)	    "
* phitci:	"	"	"     (partie imaginaire)    	    "
* dphitr : Premiere approximation de l'objet (partie reelle)      real*8	
* dphiti :	"	"	"     (partie imaginaire)	    "
* Nor_phi_T : NORM de la premiere approximation	      	        real*8
* phitr : First approximation of the object
*
	implicit none
*
	integer*4 n,pntr,i,npix(2)
	real*4	alphat,s(n),h(n),snr(n),phicr(n),phici(n)
        real*4 	phitr(n),phitcr(n),phitci(n),nor_phi_t,ratio
	real*8	dphitr(n), dphiti(n)
*
	do i=1,n
	 if (snr(i).ge.alphat .and. h(i).ne.0.) then
	  ratio = s(i)/h(i)
	  phitcr(i) = ratio * phicr(i)
	  phitci(i) = ratio * phici(i)
	 else
	  phitcr(i) = 0.
	  phitci(i) = 0.
	 end if
	 dphitr(i) = dble(phitcr(i))
	 dphiti(i) = dble(phitci(i))
	end do
*
* Calcul de la norme de phitc
*
	call normel2(n,phitcr,phitci,nor_phi_t)
*
* Calcul de phi_t
*
	npix(1) = int(sqrt(float(n)))
	npix(2) = npix(1)
	call fourier1(dphitr,dphiti,npix(1),npix(2),-1)
* Assume zero phase (real image), so I take the modulus of the inverse FFT:
	call mod_tf(n,dphitr,dphiti)
*
	call d_to_s(dphitr,phitr,n)
*
	return
	end
*
C----------------------------------------------------------------------
	subroutine dnormel1(nx,ny,a,b)
*
* Le tableau B contiendra les valeurs normalisees L1 du tableau A
*
	implicit none
*
	integer*4	nx,ny
	real*8		a(nx,ny), b(nx,ny)
*
 	integer*4	i, j
	real*8		somme
*
* Calcul de la somme des pixels
*
	somme = 0.
	do j = 1, ny
	 do i = 1, nx
	  somme = somme + a(i,j)
	 end do
	end do
*
* Normalisation
*
	do j = 1, ny
	 do i = 1, nx
	  b(i,j) = a(i,j)/somme
	 end do
	end do
*
	return
	end
*
*
	subroutine normel2(n,reel,imagi,norm)
*
* Calcul de la norme L2 du complexe (REEL, IMAGI) dans l'espace de Fourier
*
	integer*4 n, i
	real*4		reel(*), imagi(*)
	real*4		norm
*
*
	norm = 0.
	do i = 1, n
	 norm = norm+(reel(i)*reel(i)+imagi(i)*imagi(i))
	end do
*
	norm = sqrt(norm/float(n))
*
	return
	end
*
C----------------------------------------------------------------------
C	include 'd_utilities.for'
C Contains recentre:
C	include 'jlp_fft.for'
