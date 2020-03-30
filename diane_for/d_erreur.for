*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*	PROGRAM D_ERREUR
*
* Calcul de Epsilon, Epsilon_i et Epsilon_o qui seront utilises par "cgradient"
* (stored in the symbol "D_EPSILON")
*
* INPUT:
*	- Fonction de regularisation 	*.G
*	- Fonction de lissage        	*.FLI
*	- Fonction de transfert bornee  *.FTB
*	- Rapport Signal/Bruit       	*.SNR
*	- sigmai			*.SIG
*	- Stabilisateur			*.KR
*	- TF de la premiere approximation de l'Objet (R et I)  *.PCI, *.PCR
*
*	S. ROQUES - J. VIGNEAU				01 DECEMBRE 1989
*
* JLP
* Version 12-02-98
*-------------------------------------------------------------------------------
	program d_erreur
	implicit none
	integer*4  npix(2),p_g,p_s,p_h,p_snr,p_sig,p_pcr,p_pci,p_kr,
     1	long,N,lu,i,istatus
	real*4	epsilon (3),AlphaT
	character name*20,filesum*60,comments*80,buffer*90
	character*60 file_g,file_snr,file_sig,file_h,file_s,file_pcr,
     1	file_pci,file_kr
	integer*4 madrid(1)
	common /vmr/madrid
*
* Initialisation
*
	call jlp_begin
	call jlp_inquifmt
        write(6,102)
102     format('******* Program d_erreur  Version 12-02-98')

*
55	write(6,*) ' Generic name (should end with a dot) ?'
	read(5,10) name
10	format(a)
	long=index(name,'.')-1
	if(long.le.0)goto 55
 
	lu=1
	filesum=name(:long)//'.SUM'
C "APPEND" is not available for IBM...
C       open(lu,file=filesum,status='OLD',access='APPEND')
        open(lu,file=filesum,status='OLD',access='SEQUENTIAL')
        do i=1,10000
        read(lu,10,end=989) buffer
        enddo
989	write(lu,980)
980	format(/,' MODULE ERREUR',/,
     1	' -------------',/)
 
* Ouverture du fichier de la fonction de regularisation
*
	file_g=name(:long)//'.G'
	write(6,*) ' Reading the Regularisation function',file_g
	call jlp_vm_readimag(p_g,npix(1),npix(2),file_g,comments)
* Recenter image to get zero frequency at (1,1)
        call recentre(madrid(p_g),madrid(p_g),npix(1),npix(2),npix(1))
*
* Ouverture du fichier de la fonction de lissage
*
	file_s=name(:long)//'.FLI'
	write(6,*) ' Reading the smoothing function',file_s
	call jlp_vm_readimag(p_s,npix(1),npix(2),file_s,comments)
* Recenter image to get zero frequency at (1,1)
        call recentre(madrid(p_s),madrid(p_s),npix(1),npix(2),npix(1))
*
* Ouverture du fichier de la fonction de Transfert
*
C	file_h=name(:long)//'.FTB'
	write(6,*) ' Name of the limited Transfer function:'
        read(5,10) file_h
	call jlp_vm_readimag(p_h,npix(1),npix(2),file_h,comments)
* Recenter image to get zero frequency at (1,1)
        call recentre(madrid(p_h),madrid(p_h),npix(1),npix(2),npix(1))
*
* Ouverture du fichier du rapport Signal/Bruit
*
	file_snr=name(:long)//'.SNR'
	write(6,*) ' Reading the SNR file',file_snr
	call jlp_vm_readimag(p_snr,npix(1),npix(2),file_snr,comments)
* Recenter image to get zero frequency at (1,1)
        call recentre(madrid(p_snr),madrid(p_snr),npix(1),npix(2),npix(1))
*
* on recherche le descripteur "ALPHAT"
C Doesn't work yet (Version 02-05-91):
C	call jlp_get_symbol('D_ALPHAT',buffer,istatus)
C	read(buffer,*,err=87) alphat
C	if(istatus .eq. 0) goto 88
87	write(6,*) ' Threshold for the SNR (AlphaT, used in "lissage")?'
	read(5,*) alphat
*
88	write(6,*) ' Threshold for the SNR (AlphaT):',alphat
*
* Ouverture du fichier de sigmai
*
	file_sig=name(:long)//'.SIG'
	write(6,*) ' Reading the sigma file',file_sig
	call jlp_vm_readimag(p_sig,npix(1),npix(2),file_sig,comments)
* Recenter image to get zero frequency at (1,1)
        call recentre(madrid(p_sig),madrid(p_sig),npix(1),npix(2),npix(1))
*
* Ouverture du fichier de Phi T chapeau (Reelle)
*
	file_pcr=name(:long)//'.PCR'
	write(6,*)' Reading the first approx. of the object FTR',file_pcr
	call jlp_vm_readimag(p_pcr,npix(1),npix(2),file_pcr,comments)
* Recenter image to get zero frequency at (1,1)
        call recentre(madrid(p_pcr),madrid(p_pcr),npix(1),npix(2),npix(1))
*
* Ouverture du fichier de Phi T chapeau (Imaginaire)
*
	file_pci=name(:long)//'.PCI'
	write(6,*)' Reading the first approx. of the object FTI',file_pci
	call jlp_vm_readimag(p_pci,npix(1),npix(2),file_pci,comments)
* Recenter image to get zero frequency at (1,1)
        call recentre(madrid(p_pci),madrid(p_pci),npix(1),npix(2),npix(1))
*
* Ouverture du fichier du Stabilisateur
*
	file_kr=name(:long)//'.KR'
	write(6,*)' Reading the stabilizing factor :',file_kr
	call jlp_vm_readimag(p_kr,npix(1),npix(2),file_kr,comments)
* Recenter image to get zero frequency at (1,1)
        call recentre(madrid(p_kr),madrid(p_kr),npix(1),npix(2),npix(1))
*
* Main routine now: 
*
	N = npix(1) * npix(2)
	call sp_erreur(n,madrid(p_g),madrid(p_s),madrid(p_h),madrid(p_snr),
     1	madrid (p_sig),madrid (p_pcr),madrid (p_pci),
     1	madrid (p_kr),alphat, epsilon)
*
	write(buffer,*) (epsilon(i),i=1,3)
	call jlp_set_symbol('D_EPSILON',buffer,istatus)
	write(6,101) (epsilon(i),i=1,3)
	write(lu,101) (epsilon(i),i=1,3)
101	format(' epsilon (squared): eps1, eps2, eps3 : ',3(G12.4,1x))
*
* Sortie normale
	close(lu)
	call jlp_end
	stop
	end
*
*-------------------------------------------------------------------------------
*
 	subroutine sp_erreur (n, g, s, h, snr, sigmai, phitcr, phitci,
     1		kr, alphat, epsilon)
*
* EPSILON (1) = (EPSILON(2)**2 + EPSILON(3)**2)
*
*	Epsilon(2) = Epsilon_i = ||G*St*sigmai||
*	Epsilon(3) = Epsilon_o = ||kr*S*sigmao||
*
*  On prend sigmao = | PhitC |
*
	implicit none
	integer*4 n,i
	real*4	g(n),s(n),snr(n),h(n),sigmai(n),x,
     1	phitci(n),phitcr(n),kr(n),epsilon(3),alphat
*
* Calcul de epsilon_i
*
	epsilon(2)=0.
	do i=1,n
	 if (snr(i).ge.alphat .and. h(i).ne.0.) then
	  x=g(i) * (s(i)/h(i)) * sigmai(i)
	 else
	  x=0.
	 end if
	 epsilon(2)=epsilon(2) + x*x
	end do
 	epsilon(2)=epsilon(2)/float(n)
*
* Calcul de Epsilon_o
*
	epsilon(3) = 0.
	do i=1,n
	 x=sqrt(phitcr(i)*phitcr(i) + phitci(i)*phitci(i))
	 x=x*s(i)*kr(i)
	 epsilon(3)=epsilon(3) + x*x
	end do
	epsilon(3)=epsilon(3)/float(n)
*
* calcul de epsilon (au carre)
*
	epsilon(1)=epsilon(2)+epsilon(3)
*
	return
	end
