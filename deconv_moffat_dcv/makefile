##############################################################################
# Makefile for SOURCE
# to generate executable version of JLP fortran programs
# Uses jlp_make.mk
#
# JLP
# Version 04/10/2008
##############################################################################
include $(JLPSRC)/jlp_make.mk
LBFGSLIB =/d/liblnx/math/lbfgs.a

.SUFFIXES:  .cpp .c .o .exe 
.cpp.exe:
	$(CPP) $(CFLAGS) -o $(EXEC)/$*.exe $*.cpp \
	$(JLIB) $(FITSLIB) $(MATHLIB) $(LBFGSLIB) $(F77LIB) $(XLIB) -lm

.c.o:
	$(CC) -c $(CFLAGS) $*.c

.c.exe:
	$(CC) -c $(CFLAGS) $*.c
	$(CPP) -o $(EXEC)/$*.exe $*.o \
	$(JLIB) $(FITSLIB) $(MATHLIB) $(F77LIB) $(XLIB) -lm

all: dcv_deconv_1D.exe dcv_deconv_2D.exe fit_moffat.exe clean_deconv.exe 

clean_deconv.exe : clean_deconv.c 

dcv_test0.exe : dcv_test0.c 

dcv_test1.exe : dcv_test1.c 

dcv_deconv_1D.exe : dcv_deconv_1D.c dcv_1D_sub.c dcv_1D_sub.o
	$(CC) -c $(CFLAGS) $*.c
	$(CPP) -o $(EXEC)/$*.exe $*.o dcv_1D_sub.o \
	$(JLIB) $(FITSLIB) $(MATHLIB) $(F77LIB) $(XLIB) -lm

dcv_deconv_2D.exe : dcv_deconv_2D.c dcv_2D_sub.c dcv_2D_sub.o jlp_minimize.o
	$(CC) -c $(CFLAGS) $*.c
	$(CPP) -o $(EXEC)/$*.exe $*.o dcv_2D_sub.o jlp_minimize.o \
	$(JLIB) $(FITSLIB) $(MATHLIB) $(F77LIB) $(XLIB) -lm

jlp0_rdfitsio.o : ../jlpsub0/jlp0_rdfitsio.c jlp_dcv.h
	$(CC) -c $(CFLAGS) -o $*.o ../jlpsub0/$*.c

jlp0_wrfitsio.o : ../jlpsub0/jlp0_wrfitsio.c jlp_dcv.h
	$(CC) -c $(CFLAGS) -o $*.o ../jlpsub0/$*.c

jlp_syst.o : jlp_syst.c jlp_dcv.h

jlp_minimize.o : ../sourcc/jlp_minimize.c
	$(CC) -c $(CFLAGS) -o $*.o ../sourcc/$*.c

dcv_1D_sub.o : dcv_1D_sub.c jlp_dcv.h

dcv_2D_sub.o : dcv_2D_sub.c jlp_dcv.h

fit_moffat.exe : fit_moffat.cpp

clean:
	rm -f *.o
