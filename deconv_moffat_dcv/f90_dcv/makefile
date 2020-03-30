##############################################################################
# Makefile for SOURCE
# to generate executable version of JLP fortran programs
# Uses jlp_make.mk
#
# March 2002: not fully working on Vega workstation.
#
# on bonnie, see Intel compiler: /opt/intel/compiler50/ia32
#
#
# JLP
# Version 04/04/2002
##############################################################################
include $(JLPSRC)/jlp_make.mk
DCV_OBJ=
DCV_OBJ=dcv_cgrad.o
# NB: to keep the *.f, use vf90
# For virgo and F95/NAG compiler (to avoid pb with names with underscores:
F90=f95 -f77
NAGCP=/usr/local/lib/NAGWare/libf95.a
FFLAGS=-g 
# For vega (March 2001)
F90=f95
NAGCP=
F77LIB=
MATHLIB = /usr/local/lib/libnag.a $(JLPLIB)/math/libfftw.a
# EOF vega


.SUFFIXES:  .f90 .for .o .exe 

.f90.o:
	$(F90) -c $(FFLAGS) $*.f90

.f90.exe:
	$(F90) -c $(FFLAGS) $*.f90
	$(F90) $(FFLAGS) -o $(EXEC)/$*.exe $*.o \
	$(DCV_OBJ) $(JLIB) $(MIDLIB) $(MATHLIB) $(XLIB) $(F77LIB)
#	rm $*.o

.for.o:
	$(F77) -c $(FFLAGS) $*.for

.for.exe:
#	$(EXEC)/esoext1.exe -f $*.for
	$(F77) -c $(FFLAGS) $*.for
	$(F77) $(FFLAGS) -o $(EXEC)/$*.exe $*.o \
	$(DCV_OBJ) $(JLIB) $(MIDLIB) $(MATHLIB) $(XLIB) $(F77LIB)
#	rm $*.o

#.o.exe:
#	$(F77) $(FFLAGS) -o $(EXEC)/$*.exe $*.o \
#	$(JLIB) $(MATHLIB) $(MIDLIB) $(XLIB) $(F77LIB)
#	rm $*.o

all: $*.exe
$*.exe : $*.for
#$*.exe : $*.o $*.for

dcv_test0.exe : dcv_test0.for dcv_cgrad.o

dcv_test1.exe : dcv_test1.f90 dcv_cgrad.o
	$(F90) -c $(FFLAGS) $*.f90
	$(F90) $(FFLAGS) -o $(EXEC)/$*.exe $*.o \
	$(JLIB) $(MIDLIB) $(MATHLIB) $(XLIB) $(F77LIB)

# Old version with 3 modules
dcv_deconv_1DD.exe : dcv_deconv_1D.f90 dcv_1D_data_mod.f90 \
dcv_cgrad_mod.f90 dcv_1D_mod.f90
	$(F90) -c $(FFLAGS) dcv_1D_data_mod.f90 
	$(F90) -c $(FFLAGS) dcv_1D_mod.f90 
	$(F90) -c $(FFLAGS) dcv_cgrad_mod.f90 
	$(F90) -c $(FFLAGS) dcv_deconv_1D.f90
	$(F90) $(FFLAGS) -o $(EXEC)/$*.exe dcv_deconv_1D.o dcv_1D_data_mod.o\
	dcv_cgrad_mod.o dcv_1D_mod.o $(JLIB) $(MIDLIB) $(MATHLIB) $(XLIB) \
	$(F77LIB)

# New version with 2 modules
dcv_deconv_1D.exe : dcv_deconv_1D.f90 dcv_1D_data_mod.f90 \
dcv_1D_mod.f90
	$(F90) -c $(FFLAGS) dcv_1D_data_mod.f90 
	$(F90) -c $(FFLAGS) dcv_1D_mod.f90 
	$(F90) -c $(FFLAGS) dcv_deconv_1D.f90
	$(F90) $(FFLAGS) -o $(EXEC)/$*.exe dcv_deconv_1D.o dcv_1D_data_mod.o\
	dcv_1D_mod.o $(JLIB) $(MIDLIB) $(MATHLIB) $(XLIB) $(F77LIB) $(NAGCP)

dcv_deconv_2D.exe : dcv_deconv_2D.f90 dcv_2D_data_mod.f90 \
dcv_2D_mod.f90
	$(F90) -c $(FFLAGS) dcv_2D_data_mod.f90 
	$(F90) -c $(FFLAGS) dcv_2D_mod.f90 
	$(F90) -c $(FFLAGS) dcv_deconv_2D.f90
	$(F90) $(FFLAGS) -o $(EXEC)/$*.exe dcv_deconv_2D.o dcv_2D_data_mod.o\
	dcv_2D_mod.o $(JLIB) $(MIDLIB) $(MATHLIB) $(XLIB) $(F77LIB) $(NAGCP)

dcv_cgrad.o : dcv_cgrad.for

dcv_cgrad0.o : dcv_cgrad0.f90

dcv_cgrad_test.o : dcv_cgrad_test.for
