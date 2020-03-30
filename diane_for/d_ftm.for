*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*	PROGRAM FTM
*
* Normalisation L1 de la reponse impulsionnelle
* Transformee de Fourier de la reponse impulsionnelle normalisee
*
* Input: *.PSF , reponse impulsionnelle
* Output: *.FTB, module de la Fonc. de transfert bornee
*
* S. ROQUES - J. VIGNEAU				30 NOVEMBRE 1989
* JLP (Modif Midas)
* Version 02-05-91
*-------------------------------------------------------------------------------
	program d_ftm
	implicit none
	integer*4 nx,ny,PntrIn,PntrOut,Pnorm
	integer*4 dPreal, dPimag,long,N,lu
	real*4	lim
	character Name*20,comments*80,buffer*90
	character*40 FileIn,FileOut,FileSum
	integer*4 madrid(1)
	common /vmr/madrid
*
	call jlp_begin
*
* Initialisation
*
	call jlp_inquifmt
	write(6,102)
102     format('******* Program d_ftm  Version 02-07-92')

	write(6,12)
12	format('Normalization and FFT of the PSF')
*
* on recherche le nom a donner aux fichiers
*
55	write(6,*) ' Generic name for the files (end with a dot)'
	read(5,10) name
10	format(a)
	long=index(name,'.')-1
	if(long.le.0)goto 55
*
* Lower threshold for PSF support in Fourier domain
	write(6,*) ' Lower threshold of PSF in Fourier domain (0.01?)' 
	read(5,*) lim
*
* on ouvre le fichier recapitulation
*
	lu=1
	FileSum= (name(:long)//'.SUM')
C "APPEND" is not available for IBM...
C       open(lu,file=filesum,status='OLD',access='APPEND')
        open(lu,file=filesum,status='OLD',access='SEQUENTIAL')
        do i=1,10000
        read(lu,10,end=989) buffer
        enddo
989	write (lu,980)
980	format(/,' MODULE FTM',/,' ----------')
*
* Ouverture du fichier contenant la reponse impulsionnelle
*
	FileIn=name(:long)//'.PSF'
	write(6,*) ' Reading the input PSF (centred in the frame)',FileIn
	call jlp_vm_readimag(PntrIn,nx,ny,FileIn,comments)
	write (lu, '(a)') (' Input PSF : '//FileIn)
	write (lu, '(a)') (' '//comments)
*
* on normalise la reponse impulsionnelle
*
	N = nx * ny
	call jlp_getvm (Pnorm,N*4)
	call normel1(nx,ny,madrid(PntrIn),madrid(Pnorm))
*
	call jlp_getvm (PntrOut,N*4)
* On reserve de la place pour le calcul de la T.F. en D.P.
*
	call jlp_getvm (dPreal,N*8)
	call jlp_getvm (dPimag,N*8)
*
* et on passe en double precision
*
	call S_to_D (madrid(Pnorm), madrid(dPreal), N)
*
* calculons maintenant la T.F.
*
	call fourier1(madrid(dpreal),madrid(dpimag),nx,ny,1)
*
* Module de la fonction de transfert
*
	call mod_tf(N,madrid(dPreal),madrid(dPimag))
*
* repassons en REAL*4 et liberons la memoire virtuelle
*
 	call D_to_S (madrid(dPreal), madrid(Pntrout), N)
	call jlp_freevm (dPreal,N*8)
	call jlp_freevm (dPimag,N*8)
*
* On borne la fonction de transfert
*
	call borne(n,madrid(pntrout),lim)

* Recenter image to get zero frequency at (1,1):
        call recentre(madrid(pntrout),madrid(pntrout),nx,ny,nx)
*
* et on ecrit la fonction de transfert bornee dans le fichier de sortie
*
	FileOut = (name(:long)//'.FTB')
	write(6,*)' Writing the limited Transfer Function (modulus)',
     1	fileout
	comments = 'Limited Transfer Function (modulus, centred)'
	call jlp_writeimag(madrid(pntrout),nx,ny,
     1	nx,fileout,comments)
*
* Sortie normale
*
	close(lu)
	call jlp_end
	stop
	end
*
*-------------------------------------------------------------------------------
*
	subroutine borne(n,a,lim)
*
	implicit none
*
	integer*4	n
	real*4 		A(n), lim
*
	integer*4	i
*
	do i = 1, n
	 if(a(i).le.lim) then
	  a(i) = 0.
	 end if
 	end do
*
	return
	end
*
C-----------------------------------------------------------------------------
	subroutine normel1(nx,ny,a,b)
*
* Le tableau B contiendra les valeurs normalisees L1 du tableau A
*
	implicit none
*
	integer*4	nx,ny
	real*4		a(nx,ny), b(nx,ny)
*
 	integer*4	i, j
	real*4		somme
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
C-----------------------------------------------------------------------------
C	include 'd_utilities.for'
