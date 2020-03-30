##############################################################################
# Makefile for diane 
# to generate executable version of JLP fortran programs
# JLP
# Version 12-03-91
##############################################################################
include $(JLPSRC)/jlp_make.mk
OBJ=    d_utilities.o d_cgradient.o d_erreur.o \
        d_ftmm.o d_lissage.o d_regul.o d_snr.o \
        d_fstep.o w_lissage.o
EXECUT= d_cgradient.exe d_erreur.exe d_ftm.exe d_lissage.exe d_regul.exe \
        d_repons.exe d_snr.exe
EXECUTA= $(EXEC)/d_cgradient.exe $(EXEC)/d_erreur.exe \
        $(EXEC)/d_ftm.exe $(EXEC)/d_lissage.exe $(EXEC)/d_regul.exe \
        $(EXEC)/d_repons.exe $(EXEC)/d_snr.exe $(EXEC)/w_lissage.exe
#

.SUFFIXES:
.SUFFIXES:  .o .for .exe $(SUFFIXES) 

.for.exe:
	runs esoext1 -I ../jlpsub  -f $*.for
	$(F77) -c $(FFLAGS) $*.f
	rm $*.f
	$(F77) $(FFLAGS) -o $(EXEC)/$*.exe $*.o d_utilities.o fft_jlp.o \
	$(MATHLIB) $(JLIB) $(MIDLIB) $(XLIB) $(LIBC)

.for.o:
	runs esoext1 -I ../hrsa/ -f $*.for
	$(F77) -c $(FFLAGS) $*.f
	rm $*.f
	$(F77) $(FFLAGS) -o $(EXEC)/$*.exe $*.o d_utilities.o fft_jlp.o \
	$(MATHLIB) $(JLIB) $(MIDLIB) $(XLIB) $(LIBC)

all: $(EXECUTA) 
	@sleep 1

$(EXEC)/d_actual.exe: d_actual.for

$(EXEC)/d_fstep.exe: d_fstep.for

$(EXEC)/d_cgradient.exe: d_cgradient.for

$(EXEC)/d_erreur.exe: d_erreur.for

$(EXEC)/d_ftm.exe: d_ftm.for

$(EXEC)/d_lissage.exe: d_lissage.for

$(EXEC)/d_regul.exe: d_regul.for

$(EXEC)/d_repons.exe: d_repons.for

$(EXEC)/d_snr.exe: d_snr.for

$(EXEC)/w_lissage.exe: w_lissage.for

d_utilities.o: d_utilities.for

fft_jlp.o: ../fft/fft_jlp.for
	runs esoext1 -I ../fft -f ../fft/fft_jlp.for
	mv ../fft/fft_jlp.f .
	$(F77) -c $(FFLAGS) fft_jlp.f
	rm fft_jlp.f

clear: 
	rm $(OBJ)
