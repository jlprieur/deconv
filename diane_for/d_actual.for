C++***************************************************************************
!	PROGRAM ACTUAL
!
! Ce programme permet d'actualiser les noms logiques utilises dans les modules
! de reconstruction.
! la liste des noms logiques est figee :	F_SNR		<Image>.SNR
!						F_SIGMAI	<Image>.SIG
!						F_RFT		<Image>.RFT
!						F_IFT		<Image>.IFT
!						File_x		<Image>.PSF
!						F_FTB		<Image>.FTB
!						F_LIS		<Image>.FLI
!						F_HR		<Image>.HR
!						File_x		<Image>.PHIT
!						F_PHITCR	<Image>.PCR
!						F_PHITCI	<Image>.PCI
!						F_V		<Image>.V
!						F_G		<Image>.G
!						F_KR		<Image>.KR
!
! Le nom <Image> est contenu dans le symbole NOM_DE_FICHIER, ou bien est
! donne par le qualifieur /NOM=<Image>.
!
! Pour chaque nom logique, on recherche si le nom physique <Image>.<ext>
! existe - s'il existe, on ecrit la correspondance. Si le /CONFIRM est donne
! dans la commande, on doit confirmer le nom a donner aux fichiers.
!
! Le symbole FIN_DE_PROGRAMME = 'OK' permet de controler le deroulement
! du programme a l'interieur d'une procedure de commande
!
C AUTHOR: S. ROQUES - J. VIGNEAU	
C VERSION: 17 JUILLET 1989
C-------------------------------------------------------------------------------
	PROGRAM D_ACTUAL
!
	implicit none
!
 	integer*4	status,	igot, nc, n, cli$present, Nmax	/20/
	character*20	name, logical
	character*40	File
	character*80	Output
!
	logical		eVeConfirm
!
! Initialisation de eVe
!
	call eVeInit (status)
!
	call SetCharSymbol ('FIN_DE_PROGRAMME', 1, 'KO', igot, status)
!
! Recherche de Name
!
	if (cli$present ('NOM')) then
	 call ReadCharList ('NOM', Nmax, Name, n, status)
	 if (status .ne. 0) call eVeStop ('Probleme de saisie du Nom...')
	else
	 call GetCharSymbol ('NOM_DE_FICHIER', 1, name, igot, status)
	 if (status .ne. 0) call eVeStop ('Nom des fichiers ?...')
	end if
!
 	call str$trim (name, name, n)
!
	output =
	1('Assignation de noms logiques aux fichiers '//name(:n)//'.*')
	call PutText_ST (Output, status)
!
! O.K. ?
!
	if ( .not. eveconfirm()) call eVeStop
	1 ('Donnez la commande correcte...')
!
! On commence par effacer TOUS les noms logiques
!
	call lib$delete_logical ('F_SNR')
	call lib$delete_logical ('F_SIGMAI')
	call lib$delete_logical ('F_RTF')
	call lib$delete_logical ('F_IFT')
	call lib$delete_logical ('F_FTB')
	call lib$delete_logical ('F_LIS')
	call lib$delete_logical ('F_HR')
	call lib$delete_logical ('F_PHITCR')
	call lib$delete_logical ('F_PHITCI')
	call lib$delete_logical ('F_V')
	call lib$delete_logical ('F_G')
	call lib$delete_logical ('F_KR')
!
! Assignation de F_SNR - Fichier du rapport Signal/Bruit
!
 	File = name(:n)//'.SNR'
	Logical = 'F_SNR'
	call SP_Actual (File, Logical, status)
	call str$trim (File, File, nc)
!
	if (status .ne. 0) then
	 Output = (File(:nc)//' n''existe pas ! Lancer le module SNR')
	 go to 1
	else
	 call PutText_ST
	1 (File(:nc)//' a ete cree par le module SNR', status)
	end if
!
! Assignation de F_SIGMAI
!
 	File = name(:n)//'.SIG'
	Logical = 'F_SIGMAI'
	call SP_Actual (File, Logical, status)
	call str$trim (File, File, nc)
!
	if (status .ne. 0) then
	 Output = (File(:nc)//' n''existe pas ! Lancer le module SNR')
	 go to 1
	else
	 call PutText_ST
	1 (File(:nc)//' a ete cree par le module SNR', status)
	end if
!
! Assignation de F_RTF - TF de l'Image (partie reelle)
!
 	File = name(:n)//'.RTF'
	Logical = 'F_RTF'
	call SP_Actual (File, Logical, status)
	call str$trim (File, File, nc)
!
	if (status .ne. 0) then
	 Output = (File(:nc)//' n''existe pas ! Lancer le module SNR')
	 go to 1
	else
	 call PutText_ST
	1 (File(:nc)//' a ete cree par le module SNR', status)
	end if
!
! Assignation de F_IFT - TF de l'Image (partie Imaginaire)
!
 	File = name(:n)//'.IFT'
	Logical = 'F_IFT'
	call SP_Actual (File, Logical, status)
	call str$trim (File, File, nc)
!
	if (status .ne. 0) then
	 Output = (File(:nc)//' n''existe pas ! Lancer le module SNR')
	 go to 1
	else
	 call PutText_ST
	1 (File(:nc)//' a ete cree par le module SNR', status)
	end if
!
! Assignation de F_FTB - Fichier de la fonction de Transfert Bornee
!
	call PutText_ST (' ', status)
!
 	File = name(:n)//'.FTB'
	Logical = 'F_FTB'
	call SP_Actual (File, Logical, status)
	call str$trim (File, File, nc)
!
	if (status .ne. 0) then
	 Output = (File(:nc)//' n''existe pas ! Lancer le module FTM')
	 go to 1
	else
	 call PutText_ST
	1 (File(:nc)//' a ete cree par le module FTM', status)
	end if
!
! Assignation de F_LIS - Fonction de Lissage
!
	call PutText_ST (' ', status)
!
 	File = name(:n)//'.FLI'
	Logical = 'F_LIS'
	call SP_Actual (File, Logical, status)
	call str$trim (File, File, nc)
!
	if (status .ne. 0) then
	 Output = (File(:nc)//' n''existe pas ! Lancer le module LISSAGE')
	 go to 1
	else
	 call PutText_ST
	1 (File(:nc)//' a ete cree par le module LISSAGE', status)
	end if
!
! Assignation de F_HR - Ouverture synthetique
!
 	File = name(:n)//'.HR'
	Logical = 'F_HR'
	call SP_Actual (File, Logical, status)
	call str$trim (File, File, nc)
!
	if (status .ne. 0) then
	 Output = (File(:nc)//' n''existe pas ! Lancer le module LISSAGE')
	 go to 1
	else
	 call PutText_ST
	1 (File(:nc)//'  a ete cree par le module LISSAGE', status)
	end if
!
! Assignation de F_PHITCR
!
 	File = name(:n)//'.PCR'
	Logical = 'F_PHITCR'
	call SP_Actual (File, Logical, status)
	call str$trim (File, File, nc)
!
	if (status .ne. 0) then
	 Output = (File(:nc)//' n''existe pas ! Lancer le module LISSAGE')
	else
	 call PutText_ST
	1 (File(:nc)//' a ete cree par le module LISSAGE', status)
	end if
!
! Assignation de F_PHITCI
!
 	File = name(:n)//'.PCI'
	Logical = 'F_PHITCI'
	call SP_Actual (File, Logical, status)
	call str$trim (File, File, nc)
!
	if (status .ne. 0) then
	 Output = (File(:nc)//' n''existe pas ! Lancer le module LISSAGE')
	 go to 1
	else
	 call PutText_ST
	1 (File(:nc)//' a ete cree par le module LISSAGE', status)
	end if
!
! Assignation de F_V - Fichier du support de l'Objet
!
	call PutText_ST (' ', status)
!
 	File = name(:n)//'.V'
	Logical = 'F_V'
	call SP_Actual (File, Logical, status)
	call str$trim (File, File, nc)
!
	if (status .ne. 0) then
	 Output =
	1 (File(:nc)//' n''existe pas ! Lancer le module REGULARISATION')
	 go to 1
	else
	 call PutText_ST
	1 (File(:nc)//'   a ete cree par le module REGULARISATION', status)
	end if
!
! Assignation de F_G - Fonction de regularisation
!
 	File = name(:n)//'.G'
	Logical = 'F_G'
	call SP_Actual (File, Logical, status)
	call str$trim (File, File, nc)
!
	if (status .ne. 0) then
	 Output =
	1 (File(:nc)//' n''existe pas ! Lancer le module REGULARISATION')
	 go to 1
	else
	 call PutText_ST
	1 (File(:nc)//'   a ete cree par le module REGULARISATION', status)
	end if
!
! Assignation de F_KR - Fichier du Stabilisateur
!
 	File = name(:n)//'.KR'
	Logical = 'F_KR'
	call SP_Actual (File, Logical, status)
	call str$trim (File, File, nc)
!
	if (status .ne. 0) then
	 Output =
	1 (File(:nc)//' n''existe pas ! Lancer le module REGULARISATION')
	 go to 1
	else
	 call PutText_ST
	1 (File(:nc)//'  a ete cree par le module REGULARISATION', status)
	end if
!
! sortie normale
!
1      	call SetCharSymbol ('FIN_DE_PROGRAMME', 1, 'OK', igot, status)
!
	call eVeStop (Output)
	end
!
!-------------------------------------------------------------------------------
!
	SUBROUTINE SP_ACTUAL (FILE, LOGICAL, STATUS)
!
! On recherche si le fichier "FILE" existe dans le repertoire ou l'on se
! trouve.
!	Si le fichier existe, on lui assigne le nom logique "LOGICAL"
!
	implicit none
!
	integer*4	nc,
     1		ncf,
     1		status
!
	character*(*)	File, Logical
!
	logical		present
!
!
!
	call str$trim (Logical, Logical, nc)
!
	inquire (file = File, name = File, exist = present)
!
	if ( .not. present) then
	 status = 1
	else
	 call str$trim (File, File, ncf)
	 call lib$set_logical (Logical(:nc), File(:ncf))
	 status = 0
	end if
!
	return
	end
