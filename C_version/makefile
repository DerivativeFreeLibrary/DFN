# makefile for compiling dencon source code
CC=gcc
FC=gfortran
LD=gcc
RM=rm -f
CFLAGS=-Wall -ansi -pedantic -g
# -g for debugging, -O for optimization

#*******************************************************
# Dencon code

EXEC = \
       dencon

OBJS = \
       denconmain.o \
       dencon.o \
       halton.o \
       problem.o \
       xstore.o


SRCS = \
       denconmain.c \
       dencon.c \
       halton.c \
       problem.c \
       xstore.c
#*******************************************************
# Dencon header

HDRS = \
           dencon.h \
           basic.h
#*******************************************************

all:	dencon dencon.code.tar

#*******************************************************

dencon:	$(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o dencon -lm -lgfortran

dencon.code.tar:	dencon
	cp makefile dencon.code.zip
	rm -f dencon.code.zip
	cp makefile makefile.install
	zip dencon.code.zip \
                $(SRCS) $(HDRS) $(HDRS) \
                makefile makefile.install
	tar cvf dencon.code.tar \
                $(SRCS) $(HDRS) \
                makefile makefile.install

clean: 
	$(RM) *.o
	$(RM) $(EXEC)
	$(RM) fort.2
#*******************************************************
wrap_minmax.o:	wrap_minmax.f90 
	$(FC) $(FFLAGS) -c wrap_minmax.f90

denconmain.o:	denconmain.c $(HDRS)
	$(CC) $(CFLAGS) -c denconmain.c

dencon.o:	dencon.c $(HDRS)
	$(CC) $(CFLAGS) -c dencon.c

halton.o:	halton.c $(HDRS)
	$(CC) $(CFLAGS) -c halton.c

problem.o:	problem.c $(HDRS) 
	$(CC) $(CFLAGS) -c problem.c

xstore.o:	xstore.c $(HDRS) 
	$(CC) $(CFLAGS) -c xstore.c

#*******************************************************
