BIN=bayes_mtf
FC=gfortran
FFLAGS= -finit-local-zero -fbacktrace
#FFLAGS= -O0 -finit-local-zero -fbacktrace -fcheck=all
LD=$(FC)
LIBS = 
LINK = -static

$(BIN):\
bayes_mtf.o \
ps_new.o 
	$(FC) $(LINK) -o $@ \
bayes_mtf.o \
ps_new.o \
$(LIBS)

	bayes_mtf.o ps_new.o

.SUFFICES: .o .f

.f.o:
	$(FC) $(FFLAGS) -c $<

.f.mod:
	$(FC) $(FFLAGS) -c $<

clean:
	rm *o *mod
