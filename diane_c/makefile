##################################################################
# Makefile for c programs 
#
# JLP
# Version 05-01-93
##################################################################
include $(JLPSRC)/jlp_make.mk

.SUFFIXES:
.SUFFIXES: .o .c .exe $(SUFFIXES) 

OBJ= 

.c.o:
	cc -c $(CFLAGS) $*.c

.o.exe:
	cc -o $(EXEC)/$*.exe $*.o \
	$(JLIB) $(FITSLIB) $(MATHLIB) $(F77LIB) $(XLIB) -lm

.c.exe:
	cc -c $(CFLAGS) $*.c
	cc -o $(EXEC)/$*.exe $*.o \
	$(JLIB) $(FITSLIB) $(MATHLIB) $(F77LIB) $(XLIB) -lm 

all:

rbdf.exe: rbdf.c 
	cc -c $(CFLAGS) rbdf.c
	cc -o $(EXEC)/rbdf.exe rbdf.o \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

dc_gaussian_noise.exe: dc_gaussian_noise.c

dc_cgradient.exe: dc_cgradient.c

dc_erreur.exe: dc_erreur.c

dc_erreur0.exe: dc_erreur0.c

dc_ftm.exe: dc_ftm.c

dc_gradient.exe: dc_gradient.c

dc_gradient0.exe: dc_gradient0.c

dc_lissage.exe: dc_lissage.c

dc_regul.exe: dc_regul.c

dc_sigmai.exe: dc_sigmai.c

dc_snr.exe: dc_snr.c

dc_direct_snr.exe: dc_direct_snr.c

clear:
	rm -f *.o
