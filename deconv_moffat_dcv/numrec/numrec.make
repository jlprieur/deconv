############################################################################
# Makefile to create and update numrec_c.a
#
# JLP
# Version 05-04-2002 
############################################################################
include $(JLPSRC)/jlp_make.mk
DESTIN=$(JLPLIB)/jlp/math
DESTIN=.
# 
SRC= dbrent.c df1dim.c dlinmin.c f1dim.c frprmn.c linmin.c mnbrak.c \
     nrutil.c
SRC= dbrent.c mnbrak.c dcv_cgrad.c nrutil.c 
OBJ=$(SRC:.c=.o)
################################################################

.SUFFIXES:
.SUFFIXES: .o .c $(SUFFIXES)

.c.o:
	cc -c $(CFLAGS) -I../midas/incl $*.c
	ar ruv $(DESTIN)/numrec_c.a $*.o 

all: $(OBJ)
	ranlib $(DESTIN)/numrec_c.a

dbrent.o: dbrent.c

df1dim.o: df1dim.c

dlinmin.o: dlinmin.c

f1dim.o: f1dim.c

frprmn.o: frprmn.c

linmin.o: linmin.c

mnbrak.o: mnbrak.c

nrutil.o: nrutil.c

clear:
	rm -f $(OBJ)
