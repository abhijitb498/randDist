FC = gfortran 
source = nrtype.f90 nrutil.f90 spline_mod.f90 rand_dist_v3.f90 rand.f90
all : rand  


rand : $(source)
	$(FC) -o $@  $(source)

#install: 
#	@cp htg /usr/bin/ 
#	@cp hist_mod.mod /usr/lib/gcc/x86_64-linux-gnu/9/finclude/

.PHONY : clean

clean:
	@rm -f *.o *.mod htg
