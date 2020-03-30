#!/bin/csh
############################################################
# object_cerga.com
# To generate an image
# 
# JLP
# Version 14/07/01
############################################################
# Size 
echo "128,128" >! jlp_lu5.tmp 
# Filename 
echo "s5_psf" >> jlp_lu5.tmp 
# Number of objects
echo "1" >> jlp_lu5.tmp 
########## Object #2 #########################
# Analytical object
echo "1" >> jlp_lu5.tmp 
# Airy disk
echo "12" >> jlp_lu5.tmp 
# Position (center=0,0)
echo "0,0" >> jlp_lu5.tmp 
# Diameter
echo "16" >> jlp_lu5.tmp 
# Maximum
echo "1." >> jlp_lu5.tmp 
####
$EXEC/object_cerga.exe < jlp_lu5.tmp
end:
