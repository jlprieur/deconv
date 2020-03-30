*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Program REPONS
* Creation d'un fichier image de type BDF contenant une reponse impulsionnelle
*
* La valeur de SIGMA doit etre donnee en qualifieur
* du programme a l'interieur d'une procedure de commande
*
* Pas d'ecriture dans le fichier recapitulation
*
*****---------------------------------------------------------------------------
*	S. ROQUES - J. VIGNEAU				30 NOVEMBRE 1989
* JLP (Modif Midas)
* Version 10-01-90
*-------------------------------------------------------------------------------
	program d_repons
*
	implicit none
*
	integer*4 N,nx,ny,pntr
	real*4	sigma,rayon
	character name*40,comments*80	
	integer*4 madrid(1)
	common /vmr/madrid
*
* Initialisation
*
	call jlp_begin
	call jlp_inquifmt
        write(6,102)
102	format('******* Program d_repons  Version 02-05-91')
*
* saisie de nx
*
	write(6,11)
11	format(' Size of the output image (square): nx only')
	read(5,*) nx
	ny=nx
*
* saisie de sigma
*
	write(6,12)
12	format(' Sigma of the PSF :')
	read(5,*) sigma
*
* construction de la reponse
*
	N=nx*ny*4
	call jlp_getvm(pntr,N)
	call PSF(madrid(pntr),nx,ny,sigma,rayon)
*
* Ouverture du fichier Image
*
	write (comments,1002) sigma
1002	format(' Reponse Impulsionnelle, sigma = ',f4.1)
	write(6,13)
13	format(' Output image : (should be *.PSF)')
	name=' '
	call jlp_writeimag(madrid(pntr),nx,ny,nx,
     1	name,comments)
*
* ecriture du descripteur "rayon" de la PSF
*
	write(6,1001) rayon
1001	format (' Rayon de la Reponse impulsionnelle : ', F5.2)
*
	call jlp_end
	stop
	end
*
*
 	SUBROUTINE PSF (OBJ,NX,NY,SIGMA,RAYON)
*
* Simulation de la reponse impulsionnelle gaussienne
* constante sur des cercles concentriques
*
* Output:
* RAYON
*
	implicit none
*
	integer*4 nx,ny
	real*4		obj(nx,ny), sigma, rayon,sigma2
*
 	integer*4	ix, iy, k, l, ik, il
	real*4	      	arg, min, pi
	pi=3.14159265358979323846
*
	min = .0001
	rayon = 0.
	ix = nx/2 + 1
	iy = ny/2 + 1
	sigma2 = 2 * sigma * sigma
*
	do l = 1, ny
	 il = (l - iy)*(l - iy)
	 do k = 1, nx
	  ik = (k - ix)*(k - ix)
	  arg = - float(ik +il)/sigma2
	  obj (k,l) = exp(arg)
	  if (obj(k,l) .lt. min) then
	   obj(k,l) = 0.
	  else
	   rayon = rayon + 1.
	  end if
	 end do
	end do
*
	rayon = SQRT (rayon/pi)
	return
	end
