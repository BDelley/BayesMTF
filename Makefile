BIN=bayes_mtf
FC=gfortran
FFLAGS= -finit-local-zero -fbacktrace

$(BIN): bayes_mtf.o ps_new.o 
	$(FC) $(LINK) -o $@ bayes_mtf.o ps_new.o 

.SUFFICES: .o .f

.f.o:
	$(FC) $(FFLAGS) -c $<

.f.mod:
	$(FC) $(FFLAGS) -c $<

clean:
	rm *.o *.mod
