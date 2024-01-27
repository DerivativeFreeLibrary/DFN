
FC = gfortran
RM = rm -f

#FFLAGS = -O3
FFLAGS = -g


OBJS_CS = main_constr.o CS-DFN.o problem5.o halton.o sobol.o
OBJS = main_constr.o DFN_simple.o problem5.o halton.o sobol.o

OBJS_CS = main_constr.o CS-DFN.o problem.o halton.o sobol.o
OBJS = main_constr.o DFN_simple.o problem.o halton.o sobol.o

all: cs-dfn dfn_simple clean

cs-dfn:  $(OBJS_CS) 
	$(FC) -o cs-dfn $(OBJS_CS)

dfn_simple:  $(OBJS) 
	$(FC) -o dfn_simple $(OBJS)

.SUFFIXES : .f90 .o

.f90.o: $* ; $(FC) $(FFLAGS) -c $*.f90

clean: 
	$(RM) *.o
	$(RM) *.mod

cleanall:
	$(RM) *.o
	$(RM) *.mod
	$(RM) *~
	$(RM) dfn_simple
	$(RM) cs-dfn
	$(RM) *.2
	$(RM) *.1

