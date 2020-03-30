*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Program d_snrr
* to determine the SNR in the frequency space of an experimental image,
* for which the square modulus is known.
*
* INPUT:
*   Square modulus 
* OUTPUT
*
*	- tableau SNR de nx*ny elements		*.SNR
*	- tableau du majorant du bruit en frequences spatiales  *.SIG
*	- recapitulation des parametres de SNR dans le fichier  *.SUM
*
*
* JLP (modif midas...)
* Version 02-10-91
*-------------------------------------------------------------------------------
	program d_snrr
	real*4 lth
	character*40 filein,filesum,filesnr,filesig
	integer*4 nx,ny,pntrin,pntrsnr,pntrsig
	integer*4 n,long,lu
	character name*40,comments*80
	integer*4 madrid(1)
	common /vmr/madrid
*
* Initialisation
	call jlp_begin
	call jlp_inquifmt
        write(6,102)
102     format('******* Program d_snrr  Version 02-10-91')

*
* Reading square modulus
*
	write(6,*) ' Input file with the square modulus:' 
	filein=' '
	call jlp_vm_readimag(pntrin,nx,ny,filein,comments)
	call recentre(madrid(pntrin),madrid(pntrin),nx,ny,nx)
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
980	format(/,' MODULE D_SNRR',/,' ----------',/)
	write (lu,981) filein
981	format(' Square modulus : ',A)
*
	write (6,501)
501	format (' Lower threshold of the modulus ?')
	read(5,*) lth

	write (lu,500) lth 
500	format (' Lower threshold:',G12.5) 
*
	n = nx*ny
*
* on reserve de la place en REAL*4
	nn=n*4
	call jlp_getvm(pntrsnr,nn)
	call jlp_getvm(pntrsig,nn)
*
* Calcul du Rapport Signal/Bruit dans l'Espace des Frequences
*
	call rsb(madrid(pntrin),madrid(pntrsig),madrid(pntrsnr),n,lth)
*
* Output of the images:
*
	filesnr = Name(:long)//'.SNR'
	comments = 'SNR (in the Fourier domain) of:'//fileIn
	call jlp_writeimag(madrid(pntrsnr),nx,ny,nx,
     1	filesnr,comments)
*
* ouverture du fichier .SIG
*
	filesig = Name(:long)//'.SIG'
	comments = 'Sigma (in Fourier domain) of:'//fileIn
	call jlp_writeimag(madrid(pntrsig),nx,ny,nx,
     1	filesig,comments)
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
*-------------------------------------------------------------------------------
*
 	subroutine rsb(modsq,sigma_i,snr,n,lth)
*
* Calcule le rapport Signal/Bruit
*
	real*4 lth,value,modsq_min
	integer*4 n,i
	real*4	modsq(*),sigma_i(*),snr(*)
*
	modsq_min=lth*lth
	print *,' modsq_min',modsq_min
	do i = 1,n
	 value=modsq(i)
	 if(value.gt.modsq_min) then
	   snr(i) = sqrt(value)/lth
	   sigma_i(i)=1./snr(i)
	  else
	   snr(i)=0.
	   sigma_i(i)=1.E06
	  endif
	end do
*
	return
	end
C-----------------------------------------------------------------------------
C	include 'd_utilities.for'
C Contains recentre:
C	include 'fft_jlp.for'
