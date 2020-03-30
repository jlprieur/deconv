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
echo "s3_yy0" >> jlp_lu5.tmp 
# Number of objects
echo "7" >> jlp_lu5.tmp 
########## Object #1 #########################
# Analytical object
echo "1" >> jlp_lu5.tmp 
# Dimmed disk
echo "13" >> jlp_lu5.tmp 
# Position (center=0,0)
echo "0,0" >> jlp_lu5.tmp 
# Diameter
echo "100" >> jlp_lu5.tmp 
# Maximum
echo "10." >> jlp_lu5.tmp 
# Parameters for dimmed disk
echo "1.,1." >> jlp_lu5.tmp 
########## Object #2 #########################
# Analytical object
echo "1" >> jlp_lu5.tmp 
# Airy disk
echo "12" >> jlp_lu5.tmp 
# Position (center=0,0)
echo "10,10" >> jlp_lu5.tmp 
# Diameter
echo "5" >> jlp_lu5.tmp 
# Maximum
echo "20." >> jlp_lu5.tmp 
########## Object #3 #########################
# Analytical object
echo "1" >> jlp_lu5.tmp 
# Airy disk
echo "12" >> jlp_lu5.tmp 
# Position (center=0,0)
echo "5,-8" >> jlp_lu5.tmp 
# Diameter
echo "5" >> jlp_lu5.tmp 
# Maximum
echo "18." >> jlp_lu5.tmp 
########## Object #4 #########################
# Analytical object
echo "1" >> jlp_lu5.tmp 
# Airy disk
echo "12" >> jlp_lu5.tmp 
# Position (center=0,0)
echo "-8,5" >> jlp_lu5.tmp 
# Diameter
echo "5" >> jlp_lu5.tmp 
# Maximum
echo "20." >> jlp_lu5.tmp 
########## Object #5 #########################
# Analytical object
echo "1" >> jlp_lu5.tmp 
# Airy disk
echo "12" >> jlp_lu5.tmp 
# Position (center=0,0)
echo "-6,-2" >> jlp_lu5.tmp 
# Diameter
echo "5" >> jlp_lu5.tmp 
# Maximum
echo "18." >> jlp_lu5.tmp 
########## Object #6 #########################
# Analytical object
echo "1" >> jlp_lu5.tmp 
# Airy disk
echo "12" >> jlp_lu5.tmp 
# Position (center=0,0)
echo "-3,0" >> jlp_lu5.tmp 
# Diameter
echo "5" >> jlp_lu5.tmp 
# Maximum
echo "20." >> jlp_lu5.tmp 
########## Object #7 #########################
# Analytical object
echo "1" >> jlp_lu5.tmp 
# Circular ring 
echo "5" >> jlp_lu5.tmp 
# Position (center=0,0)
echo "-20,-10" >> jlp_lu5.tmp 
# External length of axes x and y 
echo "15" >> jlp_lu5.tmp 
echo "15" >> jlp_lu5.tmp 
# Internal length of axes x and y 
echo "12" >> jlp_lu5.tmp 
echo "12" >> jlp_lu5.tmp 
# Maximum
echo "15." >> jlp_lu5.tmp 
####
$EXEC/object_cerga.exe < jlp_lu5.tmp
end:
