*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*	PROGRAM D_CGRADIENT
*
* Central reconstruction with conjugate gradients
*
* INPUT:
*	*.PHIT  First approx. of the object
*	*.PCR   First approx. of the object: real part of FT
*	*.PCI   First approx. of the object: imag. part of FT
*	*.V  	Support of the object (direct space)
*	*.G	Regularizing function
*
* OUTPUT:
*	*.PHIN	Final estimation of the object
*
* Log file: *.SUM
*
*
* Pointers :
*  ---------
* p_Phit	PhiT		partie reelle			R4 -
* p_PhitCr	PhiT chapeau	partie reelle			R4 -
* p_PhitCi	  "     "	partie imaginaire		R4 -
* p_V		support de l'Objet				R4 -
* p_G		Fonction de regularisation			R4 -
* p_phin	Objet reconstruit (module)			R4 -
* p_Zn		|	Espace de travail pour les fonctions	R4   vm
* p_Rn		|	zn et rn (Cf. These page 42)		R4   vm
* p_ddnr	|	et pour la fonction dn		R8	     vm
* p_ddnI	|					R8	     vm
* p_ddntemp	|	tableau temporaire pour dn-1	R8	     vm
* p_dphitr	PhiT		partie reelle		R8	     vm
* p_dphiti	PhiT		partie imaginaire	R8	     vm
*
* From cgradient.for, S. ROQUES - J. VIGNEAU, Version 01-12-89
*
* JLP
* Version 02-06-92
*-------------------------------------------------------------------------------
	PROGRAM D_CGRADIENT
	implicit none
	integer*4 npix(2),n,n4,n8,p_phit,p_phitcr, p_phitci
	integer*4 p_v,p_g,p_zn,p_rn,p_phin,p_ddntemp,istatus
	integer*4 p_ddnr,p_ddni,p_dphitr,p_dphiti,long,lu,i
	real*4 eps(3)
	character name*20,comments*80,filesum*40,buffer*90
 	character*40 fileout,filephit,filepcr,filepci,filev,fileg
	integer*4 madrid(1)
	common /vmr/madrid
*
* Initialisation
	call jlp_begin
        write(6,112)
112     format(/,'******* Program d_cgradient  Version 12-06-92')
	call jlp_inquifmt

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
980	format(/,' Program CGRADIENT',/,
     1	' ----------------',/)
*
* Ouverture du fichier contenant la premiere approximation de l'Objet
*
	FilePhit = name(:long)//'.PHIT'
	write(6,*)' Reading the first approx. of the object:',filephit
	call jlp_vm_readimag(p_phit,npix(1),npix(2),filephit,comments)
* Recenter image to get gravity center at (1,1): 
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
*
	call S_to_D (madrid(p_phit), madrid(p_dphitr), N)
*
* Ouverture des fichiers contenant la TF de la premiere approximation de
* l'Objet (partie reelle et partie imaginaire)
*
	FilePcr = name(:long)//'.PCR'
	write(6,*)
     1  ' Reading Fourier transf. (real) of the first approx. of the object',
     1	' (.PCR) : ',filepcr
	call jlp_vm_readimag(p_phitcr,npix(1),npix(2),filepcr,comments)
* Recenter image to get zero frequency at (1,1): 
        call recentre(madrid(p_phitcr),madrid(p_phitcr),
     1   npix(1),npix(2),npix(1))
*
	FilePci = name(:long)//'.PCI'
	write(6,*)
     1  ' Reading Fourier transf. (real) of the first approx. of the object',
     1	' (.PCI) : ',filepci
	call jlp_vm_readimag(p_phitci,npix(1),npix(2),filepci,comments)
* Recenter image to get zero frequency at (1,1): 
        call recentre(madrid(p_phitci),madrid(p_phitci),
     1     npix(1),npix(2),npix(1))
*
* Ouverture du fichier contenant le support de l'objet
*
	FileV = name(:long)//'.V'
	write(6,*)' Reading the support (.V) of the object',filev
	call jlp_vm_readimag(p_v,npix(1),npix(2),filev,comments)
* Recenter image to get gravity center at (1,1): 
        call recentre(madrid(p_v),madrid(p_v),npix(1),npix(2),npix(1))
*
* Ouverture du fichier contenant la fonction de regularisation
*
	FileG = name(:long)//'.G'
	write(6,*)' Reading regul. function (.G) ',fileg
	call jlp_vm_readimag(p_g,npix(1),npix(2),fileg,comments)
* Recenter image to get zero frequency at (1,1): 
        call recentre(madrid(p_g),madrid(p_g),npix(1),npix(2),npix(1))
*
* et lecture du descripteur Epsilon
* (doesn't work in current version 02-05-91)
C	call jlp_get_symbol('D_EPSILON',buffer,istatus)
C	read(buffer,*,err=87) (eps(i),i=1,3)
C	if(istatus.eq.0) goto 88
87	write(6,*) ' Epsilon1, eps2, eps3 (from "erreur")?'
	read(5,*) eps(1),eps(2),eps(3)
*
88	write(6,*) ' Epsilon1, eps2, eps3:'
	write(6,*) eps(1),eps(2),eps(3)
*
* reservons de la place pour travailler...
*
        N4 = N*4
        N8 = N*8
	call jlp_getvm (p_rn,N4)
	call jlp_getvm (p_zn,N8)
 	call jlp_getvm (p_ddntemp,N8)
	call jlp_getvm (p_ddnr,N8)
	call jlp_getvm (p_ddni,N8)
*
* Ouverture du fichier contenant le resultat final
	call jlp_getvm(p_phin,N4)
*
	call restit(n,npix,madrid(p_phitcr),madrid(p_phitci),
     1	madrid(p_phit),madrid(p_dphitr),madrid(p_dphiti),
     1	madrid(p_v),madrid(p_g),madrid(p_rn),madrid(p_zn),
     1	madrid(p_ddnr),madrid(p_ddni),madrid(p_ddntemp),
     1	madrid(p_phin),p_dphitr,p_dphiti,
     1	p_ddnr,p_ddni,eps,lu)
*
* et on ferme les fichiers
*
	fileout = name(:long)//'.PHIN'
	comments = 'Resultat final'
* Recenter image to get gravity center in the middle: 
        call recentre(madrid(p_phin),madrid(p_phin),npix(1),npix(2),npix(1))
	call jlp_writeimag(madrid(p_phin),npix(1),npix(2),npix(1),
     1	fileout,comments)
*
* Sortie normale
	close(lu)
	call jlp_end
	stop
	end
*------------------------------------------------------------------------------
*
	subroutine restit (n,npix,phitcr,phitci,phit,
     1	dphitr, dphiti,v,g,rn,zn,ddnr,ddni,ddntemp,
     1	phin,p_dphitr, p_dphiti,p_ddnr, p_ddnI, eps, lu)
*
*             *********************************
*
*                  RECONSTRUCTION CENTRALE
*
*             *********************************
*
*                  VALEURS PROPRES DE A*A
*
*   Application de la methode des gradients conjugues
*   Calcul du polynome dont les racines convergent vers les valeurs propres
*
*             *********************************
*
	implicit none
*
       	integer*4  N,npix(2),lu,p_dphitr, p_dphiti,p_ddnr, p_ddnI
	integer*4 k,i,m,mult(35),icn(101),itest,isig,num,nomb,nrac
	real*4	V(N),G(N),eps (3),Rn(N),phin(N),PhitcR(N), PhitcI(N),
     1	Phit(N)
 	real*4	rac(35),xabs(101),xrn(101), xdn(101),
     1	xnorfin,xx, x,xnun,pas,xnfinc,deltan, ron,
     1	tete,tetap,xmum,major
	real*8	ddnr(N), ddnI(N),dphitr(N), dphiti(N),Zn(N),ddntemp(N)
	real*8	DG2,xnorrn,psdz,omegan
        integer*4 iter_max
*
* pour eviter la limitation a .29D-38
C JLP91	real*16	test,y,xnorrnp	
	real*8	test,y,xnorrnp	
	integer*4 madrid(1)
	common /vmr/madrid
*
* majorant de 1/mu
	major=10.	
* Max number of iterations:
        iter_max=20
*
***************************************************************************
*
*
* N = npix ** 2       nombre total de pixels
* (PhitcR,PhitcI)     tf de phit - parties reelle et imaginaire
* (Phit)	     module de phit
* p_dphitr et p_dphiti  pointeurs associes en dp
* V		   support de l'objet
* G		   fonction de regularisation
* eps		 premiere estimation d'erreur (epsi)
*
***************************************************************************
*
* initialisation polynomes
* xabs = abscisses  -  xrn = polynome rn  -  xdn = polynome dn
* icn = polynome cn des changements de signe
	do k = 1, 101
	 xabs (k) = (k-1) * 0.01
	 xrn (k) = 1.
	 xdn (k) = 1.
	 icn (k) = 0
	end do
*
* initialisation gradients conjugues
* fonction de depart phi zero = v phit
	do k = 1, N
*SP
	 phin(k)=V(k)*Phit(k)	
*dp
	 dphitr(k)=dble(phin(k))
	 dphiti(k)=0.
	end do
*
* calcul de r zero = v U* g2 (phitc - phi zero chap)
* calculs de d zero et de ||rn||**2 pour n=0
*
* p_dphitr : pointeur associe a dphitr
* p_dphiti : pointeur associe a dphiti
*
	call fourier1(madrid(p_dphitr),madrid(p_dphiti),npix(1),npix(2),1)
C JLP92 
C	print *,' phitr(4),phiti(4)',dphitr(4),dphiti(4)
*
	do k = 1, N
	 DG2 = dble(G(k)) * dble(G(k))
	 dphitr(k) = (dble(phitcr(k)) - dphitr(k)) * dg2
	 dphiti(k) = (dble(phitci(k)) - dphiti(k)) * dg2
	end do
*
	call fourier1(madrid(p_dphitr),madrid(p_dphiti),npix(1),npix(2),-1)
C JLP92
C	print *,' phitr(4),phiti(4)',dphitr(4),dphiti(4)
*
	call mod_tf(n,madrid(p_dphitr),madrid(p_dphiti))
*
	xnorrn = 0.
	do k = 1, N
	 dphiti(k) = 0.
	 rn(k) = sngl(dphitr(k)) * V(k)
	 ddnr(k) = dble(rn(k))
	 xnorrn = xnorrn + ddnr(k) * ddnr(k)
	end do
*
	write(6,105)
	write(lu,105)
105	format('  N         Test',/)
*
*******************************************************************************
*
*			I T E R A T I O N
*
*******************************************************************************
*
	I = 0
600	I = I + 1
*
* on garde l'ancien ddnr dans ddntemp
	do k = 1, N
 	 ddntemp(k) = ddnr(k)
	 ddnI(k) = 0.
	end do
*
* calcul de dn chapeau
	call fourier1(madrid(p_ddnr),madrid(p_ddni),npix(1),npix(2),1)
*
* calcul de zn et de omega n
*
	do k = 1, N
	 DG2 = dble(G(k)) * dble(G(k))
	 ddnr(k) = ddnr(k) * DG2
	 ddnI(k) = ddnI(k) * DG2
	end do
*
	call fourier1(madrid(p_ddnr),madrid(p_ddni),npix(1),npix(2),-1)
* Scalar product: (d_n | A*A d_n) 
	psdz = 0.
*
	do k = 1, N
	 Zn(k) = ddnr(k) * V(k)
	 psdz = psdz + (ddntemp(k)* 1.e10) * (Zn(k) * 1.e10)
	end do
*
	if (psdz .ne. 0.) then
	 omegan = xnorrn * 1.e20/psdz
	else
	 write(6,109)
109	 format(' d_cgradient/Fatal error when computing Omega_n'/,
     1   ' Emergency exit *** (d_n | A*A d_n) = 0')
	 write(lu,103)
	 stop
	end if
*
* calculs de phi n+1 et de ||phi n+1||**2
* calculs de r n+1 et de ||r n+1||**2
*
	xnorfin = 0.
	xnorrnp = 0.
*
	do k = 1, N
	 phin(k) = phin(k) + omegan * ddntemp(k)
	 x = phin(k) * phin(k)
	 xnorfin = xnorfin + x
	 Rn(k) = Rn(k) - sngl(omegan * Zn(k))
C JLP91	 y = qext(Rn(k)) * qext(Rn(k))
	 y = dble(Rn(k)) * dble(Rn(k))
	 xnorrnp = xnorrnp + y
	end do
*
* formation des polynomes Rn et Cn
* nombre de degres de liberte : icn(95)
* isig correspond au polynome Sn(X)
*
	itest = 0
	do k = 1, 101
	 isig = sign (1.,xrn(k))
	 xx = sngl(omegan) * xabs(k) * xdn(k)
	 xrn(k) = xrn(k) - xx
*
* ... on evite le point d'accumulation en 1 :
*
	 if (k .le. 95) then
	  if(isig .ne. sign(1.,xrn(k))) then
	   itest = 1
	   icn(k) = icn(k) + 1
	  end if
	 end if
	end do
*
* test d'arret
*
	test = xnorrnp/xnorfin
	test = major * sqrt(test)
*
	write (6,102) i, test
	write (lu,102) i, test
102	format(1x,i2,5x,g12.6)
*
	if (test.le.1.e-03.or.i.ge.iter_max) then
*
* ... itest=0 signifie que Rn n'a change de signe en aucun point xk
*
	 if (itest.eq.0.or.i.ge.iter_max) goto 500
	end if
*
* calcul de nu n
	xnun = sngl(xnorrnp/xnorrn)
*
* calcul de d n+1
	do k = 1, N
	 ddnr(k) = dble(Rn(k)) + dble(xnun) * ddntemp(k)
	end do
*
* calcul du polynome D n+1
	do k = 1, 101
	 xdn(k) = xrn(k) + xnun * xdn(k)
	end do
*
* reinitialisations
	xnorrn = dble(xnorrnp)
*
* retour a I = I + 1
*
	if(i.lt.iter_max) goto 600
*
* ********************  FIN D'ITERATION   ********************
*
500	continue
*
	write(6,1000) icn(95)
	write(lu,1000) icn(95)
1000	format (/,' Number of degrees of freedom : ', i3,/)
*
* calcul des valeurs propres de A*A
*
* rac est la racine de rn
* num est son numero d'ordre
* mult est sa multiplicite
* nomb est le nombre de racines (y compris la multiplicite)
*
	pas = 1./(2.*100.)
	m = 94
	num = 0
	nomb = 0
*
	write(6,106)
	write(lu,106)
106	format('  N   Eigenvalues       Multiplicity',/)
	do k = 1, m
	 nrac = icn(k+1) - icn(k)
	 if (nrac .ne. 0) then
	  num = num + 1
	  rac (num) = xabs(k+1) - pas
	  mult (num) = nrac
	  nomb = nomb + nrac
	  write (lu,103) num,rac(num),mult(num)
	  write (6,103) num,rac(num),mult(num)
103	  format(1x, i2, 10x, f6.4, 10x, i2)
	 end if
	end do
*
* calcul des majorants de l'erreur
*
	do k = 1, N
	 dphitr(k) = dble(G(k)) * (dphitr(k) - phin(k))
	end do
*
	call dnormel2 (N, dphitr, dphiti, deltan)
	xnfinc = 0.
	do k = 1, N
	 xnfinc = xnfinc + phin(k) * phin(k)
	end do
	xnfinc = sqrt (xnfinc)
	ron = sqrt( abs(eps(1) - deltan*deltan))
*
* borne superieure de l'erreur
*
	tete = sqrt(rac(1) * xnfinc)
	tete = ron / tete
*
* calcul du 1/mu moyen
* calcul des ereurs correspondantes
*
	xmum = 0.
	do k = 1, num
	 xmum = xmum + (1./rac(k))*mult(k)
	end do
*
	xmum = xmum/float(nomb)
	xmum = (xmum + 1.)/2.
	xmum = sqrt(xmum)
	tetap = (xmum*ron)/xnfinc
*
	write (6, 1001) xnfinc,xmum,sqrt(eps(1))
	write (lu, 1001) xnfinc,xmum,sqrt(eps(1))
1001	format(/,' Norm of the solution Phi n  ( ||Phi n|| ) :',
     1	'. . . . ',g12.5,/,
     1	' Weighted mean of the inverses of the eigenvalues (1/mu): ',
     1	g12.5,/,' Error of image reconstruction (sqrt(eps1)) :',
     1	' . . . . ',g12.5)
*
	write (6, 1002) sqrt(eps(3)),sqrt(eps(2)),deltan,tetap
	write (lu, 1002) sqrt(eps(3)),sqrt(eps(2)),deltan,tetap
1002	format(' Error due to the choice of Hr: ',
     1	'(epsilon_o) sqrt(eps3) : ',g12.5,/,
	1 	' Erreur related to the SNR ratio: ',
     1	'(epsilon_i) sqrt(eps2): ',g12.5,/,
     1	' Term to correct the error (Delta_n) :',
     1	' . . . .  ',g12.5,/,
	1 	' Relative error for object reconstruction:',
     1	'  . . . .  ',g12.5)
*
	return
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
