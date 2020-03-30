*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*	PROGRAM D_REGUL(ARISATION)
*
* - Prelevement du support de l'Objet,
* - calcul de sa mesure (tau),
* - calcul de la fonction de regularisation g,
* - calcul de la mesure de 1 - g**2 (nu),
* - calcul du parametre d'interpolation (eta),
* - test sur le bon conditionnement du probleme.
*
* INPUT:
*	*.PHIT  First approx. of the object
*	*.snr	SNR of the input image (direct space)
*	*.hr	Synthetic aperture
* OUTPUT:
*	*.V  	Support of the object (direct space)
*	*.G	Regularizing function
*	*.KR	Stabilizing function
*
* With an isophote level for defining the support: 
* runs w_regul alpha_prime_T generic_name isophote_value alpha_T 
*
* With a mask contained in a file:
* runs w_regul alpha_prime_T mask.fits value alpha_T 
*
* S. ROQUES - J. VIGNEAU	14 SEPTEMBRE 1989
*
* 11-01-90: Modif. Midas
* JLP
* Version 02-05-91
*-------------------------------------------------------------------------------
	program d_regul
	implicit none
	integer*4 p_phit,p_v,p_snr,p_hr,p_g,p_kr,p_mask,long
	integer*4 lu,n,istatus,mask_in_file
	integer*4 nx,ny,nx1,ny1
	real*4 niveau,tau,alphapt,alphat
	character comments*80,name*20,buffer*90
	character*40 filephit,filev,filesnr,fileg,filekr,filehr,filesum
	character mask_name*40
	integer*4 madrid(1)
	common /vmr/madrid
*
* Initialisation
	call jlp_begin
	call jlp_inquifmt
        write(6,112)
112     format('******* Program d_regul  Version 02-05-91')

*
* Saisie du qualifieur ALPHA PRIME T (ALPHAPT)
	write(6,*) ' Threshold for the reliability of spectral',
     1	' information alphapt=4 ?'
	read(5,*) alphapt
*
55	write(6,*) ' Generic name (should end with a dot) ?'
	read(5,10) name
10	format(a)
	long=index(name,'.')-1
	if(long.le.0)goto 55
*
* ouverture du fichier de la recapitulation
	lu=1
	filesum=name(:long)//'.SUM'
C "APPEND" is not available for IBM...
C       open(lu,file=filesum,status='OLD',access='APPEND')
        open(lu,file=filesum,status='UNKNOWN',access='SEQUENTIAL')
        do i=1,10000
        read(lu,10,end=989) buffer
        enddo
989	write (lu,980)
980	format(/,' MODULE REGULARISATION',/,
     1	' ---------------------')
*
* Saisie du parametre Niveau
        mask_in_file = 1	
	write(6,56) 
56	format(' Threshold of the isophote for the support ?',
     1  '(or mask file name)')
	read(5,10) buffer 
        niveau=0.
	read(buffer,*,err=57) niveau
C IBM assumes 0. if error, and does not go to label 57
        if(niveau.eq.0.)goto 57
        mask_in_file = 0	
	write (6, 100) niveau
	write (lu, 100) niveau
100	format (' Threshold of the support : ',F9.4)
57        if(mask_in_file .eq. 1) then
	  read(buffer,10) mask_name
* Opening the mask file:
*
	  write(6,*) ' Reading the mask file ',mask_name
	  call jlp_vm_readimag(p_mask,nx,ny,mask_name,comments)
	  write (lu, 110) mask_name 
110	format (' Input mask : ',A40)
C Recentering the image to get zero frequency at (0,0):
          call recentre(madrid(p_mask),madrid(p_mask),nx,ny,nx)
	endif
*
* Ouverture du fichier contenant phit:
*
	filephit = name(:long)//'.PHIT'
	write(6,*) ' Reading the first approx. of the object',filephit
	call jlp_vm_readimag(p_phit,nx1,ny1,filephit,comments)
C Recentering the image to get zero frequency at (0,0):
        call recentre(madrid(p_phit),madrid(p_phit),nx1,ny1,nx1)
	if(mask_in_file.eq.1)then
          if(nx.ne.nx1.or.ny.ne.ny1)then
	    print *,' d_regul/ Fatal error'
            print *,' Incompatible size between .PHIT and mask!'
	    stop
	  endif
	else
	  nx=nx1
	  ny=ny1
	endif
*
* Ouverture du fichier V contenant le support de l'objet
*
	N=nx*ny
	call jlp_getvm(p_v,N*4)
*
* Fabrication du support
	if(mask_in_file .eq. 1)then
	  call sp_support1(madrid(p_mask),nx,ny,madrid(p_v),tau)
	else
	  call sp_support(nx,ny,niveau,madrid(p_phit),madrid(p_v),tau)
	endif
*
* Ouverture du fichier V contenant le support de l'objet
*
	fileV = name(:long)//'.V'
	write(6,*) ' Writing the support ',fileV
	comments = 'Support (direct space)'
C Recentering the image to get zero frequency in the center:
        call recentre(madrid(p_v),madrid(p_v),nx,ny,nx)
	call jlp_writeimag(madrid(p_v),nx,ny,nx,
     1	fileV,comments)
*
* on ouvre le fichier *.snr
*
	filesnr = name(:long)//'.SNR'
	call jlp_vm_readimag(p_snr,nx1,ny1,fileSNR,comments)
C Recentering the image to get zero frequency at (0,0):
        call recentre(madrid(p_snr),madrid(p_snr),nx1,ny1,nx1)
        if(nx.ne.nx1.or.ny.ne.ny1)then
	  print *,' d_regul/ Fatal error'
          print *,' Incompatible size between .PHIT and .SNR !'
	  stop
	endif
 
* on lit le descripteur ALPHAT dans snr
C Doesn't work in current version (02-05-91)
C	call jlp_get_symbol('D_ALPHAT',buffer,istatus)
C	print *,' istatus',istatus
C	read(buffer,*,err=87)alphat
C        if (istatus .eq. 0) goto 88
87	write(6,*) ' Threshold for the snr (alphat, used in "lissage")?'
	read(5,*) alphat
*
88	write(6,*) ' Threshold for the snr (alphat):',alphat
*
* on ouvre le fichier hr
*
	filehr = name(:long)//'.HR'
	call jlp_vm_readimag(p_hr,nx1,ny1,filehr,comments)
C Recentering the image to get zero frequency at (0,0):
        call recentre(madrid(p_hr),madrid(p_hr),nx1,ny1,nx1)
*
* Ouverture du fichier de la fonction de regularisation
	call jlp_getvm(p_g,N*4)
*
* Ouverture du fichier du Stabilisateur
	call jlp_getvm(p_kr,N*4)
*
* calcul de g, de nu et de kr
*
	call sp_g(nx,ny,madrid(p_snr),madrid(p_hr),
     1        madrid(p_g),madrid(p_kr),alphat,alphapt,tau,lu)
*
* Ouverture du fichier de la fonction de regularisation
	fileG = name(:long)//'.G'
	comments = 'Fonction de regularisation'
C Recentering the image to get zero frequency in the center:
        call recentre(madrid(p_g),madrid(p_g),nx,ny,nx)
	call jlp_writeimag(madrid(p_g),nx,ny,nx,
     1	fileG,comments)
*
* Ouverture du fichier du Stabilisateur
	filekr = name(:long)//'.KR'
	comments = 'Stabilisateur'
C Recentering the image to get zero frequency in the center:
        call recentre(madrid(p_kr),madrid(p_kr),nx,ny,nx)
	call jlp_writeimag(madrid(p_kr),nx,ny,nx,
     1	filekr,comments)
*
* Sortie normale
	close(lu)
	call jlp_end
        stop 
	end
*
*-------------------------------------------------------------------------------
*
	subroutine sp_support(nx,ny,min,phit,v,tau)
*
	implicit none
	integer*4 nx,ny
	real*4	phit(nx,ny),v(nx,ny)
	real*4  tau,min
	integer*4 i,j
*
* calcul de la mesure de V
*
	tau = 0.
*
	do j = 1, ny
	 do i = 1, nx
	  v (i, j) = 0.
	  if (phit(i, j) .ge. min) then
	   v(i, j) = 1.
	   tau = tau + 1.
	  end if
	 end do
	end do
*
	tau = sqrt (tau)
*
	return
	end
*************************************************************
* New version of sp_support with a mask in a file:
* The mask should be set to 1 on the support and to 0 outside...
* JLP91
*
*************************************************************
	subroutine sp_support1(mask,nx,ny,v,tau)
*
	implicit none
	integer*4 nx,ny
	real*4	mask(nx,ny),v(nx,ny)
	real*4  tau
	integer*4 i,j
*
* calcul de la mesure de V
*
	tau = 0.
*
        print *,' sp_support1/ nx, ny =',nx,ny
	do j = 1, ny
	 do i = 1, nx
	   v(i,j) = mask(i,j) 
	   tau = tau + v(i,j) 
	 end do
	end do
*
	tau = sqrt (tau)
*
	return
	end
*
*-------------------------------------------------------------------------------
*
	subroutine sp_g (nx,ny,snr,hr,g,kr,alphat,alphapt,tau,lu)
*
	implicit none
	integer*4 i, j, N, NM,nx,ny,lu
	real*4	snr(nx,*),hr(nx,*),G(nx,*),kr(nx,*)
        real*4 alphat,alphapt,tau,du,X,Y,Z,nu,eta
*
	N = nx
	NM = N - 1
*
* calcul de la fonction de regularisation g et du stabiliseur kr
*
	do j = 1, ny
	 do i = 1, nx
	  kr (i,j) = 0.
	  g (i,j) = 1. - hr (i,j)
	  if (G(i,j) .eq. 0.) then
	   if (snr (i,j) .ge. alphapt) then
	    G (i,j) = 1.
	   else
	    if (snr (i,j) .ge. alphat)
	1   G (i,j) = (snr (i,j) - alphat)/(alphapt - alphat)
	   end if
	  else
	   if (snr (i,j) .lt. alphat) kr (i,j) = 1.
	  end if
	  G (i,j) = 1. - G(i,j)*G(i,j)
	 end do
	end do
*
* calcul de la mesure du support d'interpolation frequentiel (1-G**2)
*
       	du = 1./float(nx)
	X = G (1,1) + G (1,N) + G (N,N) + G (N,1)
	X = X * du * du /4.
	Y = 0.
*
	do i = 2, NM
	 Y = Y + G (i,1) + G (i,N)
	end do
*
	do i = 2, NM
	 Y = Y + G (1,i) + G (N,i)
	end do
*
	Y = Y * du * du / 2.
*
	Z = 0.
*
	do j = 2, NM
	 do i = 2, NM
	  Z = Z + G(i,j)
	 end do
	end do
*
	Z = Z * du * du
*
	nu = sqrt (X + Y + Z)
	eta = tau * nu
*
* on restitue G
*
	do j = 1, ny
	 do i = 1, nx
	  G(i,j) = sqrt (abs(1 - G(i,j)))
	 end do
	end do
*
* Messages pour l'utilisateur
*
	write(6,100) tau
	write(lu,100) tau
100	format(' Mesure du support de l''Objet :          ',f5.2)
*
	write(6,101) nu
	write(lu,101) nu
101	format(' Mesure de la fonction d''interpolation : ',f5.2)
*
	write(6,102) eta
	write(lu,102) eta
102	format (' Parametre d''interpolation :             ',f5.2)
*
* Test sur le parametre d'interpolation
*
C JLP91
C	if (eta .le. 3.4) then
	if (eta .le. 6.0) then
	  write(6,103)
	  write(lu,103)
103	 format(' Le Probleme est bien conditionne',/,
     1	' on peut engager la reconstruction iterative')
	else
C	 if (eta .lt. 4.2) then
	 if (eta .lt. 8.5) then
	  write(6,104)
	  write(lu,104)
104	format(' Le probleme est moyennement stable',/,
     1	' L''erreur sur la reconstruction risque d''etre importante'
     1	,/,' on peut engager la reconstruction iterative ou bien',/,
     1	'on peut reduire le support de l''objet ou',
     1	' le gain en resolution')
	 else
	  write(6,105)
	  write(lu,105)
105	format(' Le Probleme est mal pose - la solution est instable',/,
     1	' Il faut reduire le support de l''objet ou',
     1	' le gain en resolution')
	 end if
	end if
*
	return
*
	end
