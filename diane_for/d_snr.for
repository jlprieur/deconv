*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Program SNR
* to determine the SNR in the frequency space of an experimental image,
* for which the SNR is known (in direct space).
*
* INPUT:
*   Image experimentale
*   snr0 (scalaire) : rapport Signal/Bruit dans l'espace reel
*		 (signal moyen / variance)
*   alpha (scalaire) : fiabilite de snr0 
*				alpha =	1 si snr0 est tres bien determine
*					2 valeur par defaut
*					3 si snr0 est tres mal determine.
* OUTPUT
*
*	- tableau SNR de nx*ny elements		*.SNR
*	- tableau du majorant du bruit en frequences spatiales  *.SIG
*	- TF de l'image experimentale	(partie reelle)		*.RFT
*	-	"	"	"	(partie imaginaire)	*.IFT
*
*	- recapitulation des parametres de SNR dans le fichier  *.SUM
*
*
* S. ROQUES - J. VIGNEAU				16 NOVEMRE 1989
* JLP (modif midas...)
* Version 02-05-91
*-------------------------------------------------------------------------------
	program d_snr
	implicit none
        real*4 snr0
	character*40 filein,filesum,filesnr,filesig, iftmre, iftmim
	integer*4 nx,ny,pntrin,pntrsnr,pntrsig,pntir, pntii
	integer*4 dpin,dptir,dptii,alpha,n,long,lu
	character name*20,comments*80
	integer*4 madrid(1)
	common /vmr/madrid
*
* Initialisation
	call jlp_begin
	call jlp_inquifmt
        write(6,102)
102     format('******* Program d_snr  Version 02-07-92')

*
* Saisie de snr0 
* snr0 (scalaire) : rapport Signal/Bruit dans l'espace reel
*		 (signal moyen / variance)
	write(6,11)
11	format(' Mean SNR in the direct space ?')
	read(5,*) snr0
*
* Estimation de la fiabilite
	write(6,12)
12	format(' Reliability factor of SNR (1=good, 2=medium, 3=bad) ?')
	read(5,*) alpha
*
* Ouverture du fichier Image Experimentale
*
	write(6,*) ' Experimental image ?'
	filein=' '
	call jlp_vm_readimag(pntrin,nx,ny,filein,comments)
*
* Fabrication de la chaine de caracteres "NAME"
*
55	write(6,*) ' Generic name for the output ? (dot at the end)'
	read(5,10) name
	long=index(name,'.')-1
	if(long.le.0)goto 55
10	format(a)
*
* Ouverture du fichier <name>.SUM
*
	lu=1
	filesum = Name(:long)//'.SUM'
	open(1,file=filesum,status='unknown',access='sequential')
*
* Prenons quelques notes...
*
	write (lu,980)
980	format(/,' MODULE SNR',/,' ----------',/)
	write (lu,981) filein
981	format(' Image experimentale : ',A)
*
	write (lu,500) snr0,alpha
	write (6,500) snr0,alpha
500	format (' Rapport Signal/Bruit : ',e12.5,' criterium : ',i2)
*
	N = nx*ny
*
* on reserve de la place en REAL*4
	call jlp_getvm (PntrSnr,N*4)
	call jlp_getvm (PntrSig,N*4)
  	call jlp_getvm (PnTIR,N*4)
* pour partie reelle TF image
	call jlp_getvm (PnTII,N*4)
* pour partie imaginaire TF image
 
* on reserve de la place en REAL*8
* pour l'image experimentale
	call jlp_getvm (DpIn,N*8)
* pour partie reelle TF image
  	call jlp_getvm (DpTIR,N*8)
* pour partie imaginaire TF image
	call jlp_getvm (DpTII,N*8)
*
* passage en D.P.
	call S_to_D (madrid(PntrIn), madrid(DpIn), N)
 
* on charge l'image experimentale dans DpTIR
	call psi_to_phi(nx,ny,madrid(DpIn),madrid(DpTIR))
*
* et on calcule la T.F. de l'image experimentale
*
	call fourier1(madrid(dptir),madrid(dptii),nx,ny,1)
*
* on repasse en R*4 pour la TF de l'Image Experimentale
*
	call D_to_S (madrid(DpTIR), madrid(PnTIR), N)
	call D_to_S (madrid(DpTII), madrid(PnTII), N)
*
* appel de BRUIT
* computes the sigma in PntrSig
	call bruit(nx,ny,snr0,alpha,madrid(DpIn),madrid(PntrSig))
*
* Calcul du Rapport Signal/Bruit dans l'Espace des Frequences
*
	call rsb(madrid(PntrSig),madrid(PnTIR),
     1	madrid(PnTII),madrid(PntrSnr),N)
*
* Output of the images:
*
* Recenter image to get zero frequency at (1,1):
        call recentre(madrid(pntrsnr),madrid(pntrsnr),nx,ny,nx)
	fileSnr = Name(:long)//'.SNR'
	comments = 'SNR (in the Fourier domain) of:'//fileIn
	call jlp_writeimag(madrid(pntrsnr),nx,ny,nx,
     1	filesnr,comments)
*
* ouverture du fichier .SIG
*
* Recenter image to get zero frequency at (1,1):
        call recentre(madrid(pntrsig),madrid(pntrsig),nx,ny,nx)
	filesig = Name(:long)//'.SIG'
	comments = 'Sigma (in Fourier domain) of:'//fileIn
	call jlp_writeimag(madrid(pntrsig),nx,ny,nx,
     1	filesig,comments)
*
* ouverture des fichiers .RFT et .IFT (T.F. Image Experimentale)
*
* Recenter image to get zero frequency at (1,1):
        call recentre(madrid(pntir),madrid(pntir),nx,ny,nx)
	iftmre = Name(:long)//'.RFT'
	comments = ('F.T (real part) of: '//fileIn)
	call jlp_writeimag(madrid(pntir),nx,ny,nx,
     1	iftmre,comments)
*
* Recenter image to get zero frequency at (1,1):
        call recentre(madrid(pntii),madrid(pntii),nx,ny,nx)
	iftmim = Name(:long)//'.IFT'
	comments = ('F.T (imag. part) of: '//fileIn)
	call jlp_writeimag(madrid(pntii),nx,ny,nx,
     1	iftmim,comments)
*
	close(lu)
*
*
* Sortie normale
*
	call jlp_end
	stop
 	end
*
********************************************************************
* Subroutine BRUIT
*
* To determine the noise in the Fourier domain
* The trick is to solve in 5 iterations the equation:
*      PHI_experimental(i,j) = PHI_ideal(i,j) + noise(i,j)
* The noise is assumed to be consistent with the SNR value give by the user
* and following a gaussian law.
*
* Then the Fourier transform FT(noise)(u,v) is computed from the
* map noise(i,j),
*
* INPUT: NPIX, N, K, ALPHA, PHI
* OUTPUT: SIGMAI (noise of the input image)
*
* Note:  PHI is also modified in the loop
*
********************************************************************
	subroutine bruit (nx,ny,snr0,alpha,phi,sigmai)
	implicit none
	integer*4 nx,ny,N,alpha
	real*8	phi(nx,*),alk
	real*4	SigmaI(nx,*),snr0
	integer*4 DpBr, DpBi, dpsi, i, j
	integer*4 madrid(1)
	common /vmr/madrid
*
* Setting to zero SigmaI:
	do j=1,ny
	  do i=1,nx
	    SigmaI(i,j)=0.
	  end do
	end do
 
* on reserve de la memoire
*
	N = nx*ny
	alk = float(alpha)/snr0
* Tableau R*8 du Bruit (reel)
	call jlp_getvm (dpBr,N*8)
* Tableau R*8 du Bruit (imaginaire)
	call jlp_getvm (dpBi,N*8)
* working space for the noised image
	call jlp_getvm (dpsi,N*8)
*
* initialisation du generateur gaussien (NAG)
C	call G05CBF(0)
	call jlp_random_init(1) 
*
* BOUCLE DE TRAVAIL
	do i = 1,5
*
* Computing the Gaussian noise:
* Br: gaussian noise corresponding to be present value of phi
* Bi is set to zero
	 call B_gauss (nx,ny,alk,phi,madrid(dpBr),madrid(dpBi))
*
* et on l'ajoute au signal experimental: psi=phi+Br
	 call add_bruit(nx,ny,phi,madrid(dpsi),madrid(dpBr))
*
* on calcule la TF du Bruit B
	 call fourier1(madrid(dpbr),madrid(dpbi),nx,ny,1)
*
* calcul de sigmai
* SigmaI=SigmaI+sqrt(Br**2+Bi**2)
	 call somme_i(nx,ny,sigmai,madrid(dpbr),madrid(dpbi))
*
* on remplace le signal initial par le signal bruite boucle phi = psi
	 call psi_to_phi(nx,ny,madrid(dpsi),phi)
*
	end do
*
* calcul definitif de SigmaI:
	do j = 1, ny
	 do i = 1, nx
	  sigmai(i,j) = sigmai(i,j) * 2./5.
	 end do
	end do
*
	return
	end
*
*-------------------------------------------------------------------------------
*
	subroutine b_gauss(nx,ny,alk,phi,br,bi)
*
* calcule le bruit gaussien B(i,j), real and imaginary
*
*IN: nx,ny,alk, phi
*OUT: Br, Bi
*
        real*4 work
	integer*4 nx,ny,i,j
	real*8	phi(nx,*),Sig_B,Br(nx,*), Bi(nx,*),alk
*
* calcul du bruit gaussien
*
	do j = 1,ny
	 do i = 1,nx
	  sig_B = phi(i,j)*alk
C NAG routine for random generator (G05DDF):
C	  Br(i,j) = G05DDF(phi(i,j),sig_B) - phi(i,j)
C jlp_random_gauss has a mean of 0. and a sigma of 1.
          call jlp_random_gauss(work)
	  Br(i,j) = work*sig_B
	  Bi(i,j) = 0.
	 end do
	end do
*
	return
	end
*
*-------------------------------------------------------------------------------
*
 	subroutine add_bruit (nx,ny,phi,psi,bruit)
*
	implicit none
	integer*4	nx,ny,i,j
	real*8	phi(nx,*),psi(nx,*),bruit(nx,*)
*
	do j = 1,ny
	 do i = 1,nx
	  Psi(i,j) = Phi(i,j) + Bruit(i,j)
	 end do
	end do
*
	return
	end
*
*-------------------------------------------------------------------------------
*
	subroutine somme_i (nx,ny,sigmai,b_real,b_imag)
*
	integer*4 nx,ny,i,j
	real*4	sigmai(nx,*),u,v
	real*8	b_real(nx,*),B_imag (nx,*)
*
	do j = 1, ny
	 do i = 1, nx
	  u = sngl (B_real(i,j) * B_real(i,j))
	  v = sngl (B_imag(i,j) * B_imag(i,j))
	  SigmaI(i,j) = SigmaI(i,j) + sqrt (u + v)
	 end do
	end do
*
	return
	end
*
*-------------------------------------------------------------------------------
*
	subroutine psi_to_phi(nx,ny,psi,phi)
*
	implicit none
	integer*4 nx,ny,i,j
	real*8	psi(nx,*),phi(nx,*)
*
	do j = 1, ny
	 do i = 1, nx
	  phi(i,j) = psi(i,j)
	 end do
	end do
*
	return
	end
*
*-------------------------------------------------------------------------------
*
 	subroutine rsb(sigma_i,rtf,ift,snr,n)
*
* Calcule le rapport Signal/Bruit
*
	implicit none
	integer*4 n,i
	real*4	sigma_i(n),rtf(n),ift(n),snr(n),module
*
	do i = 1,n
	 module = sqrt (rtf(i)*rtf(i) + ift(i)*ift(i))
	 snr(i) = module/sigma_i(i)
	end do
*
	return
	end
C-----------------------------------------------------------------------------
C	include 'd_utilities.for'
