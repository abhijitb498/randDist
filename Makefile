FC = gfortran 
source = nrtype.f90 nrutil.f90 spline_mod.f90 rand_mod.f90 test_rand.f90

all : rand  

rand : $(source)
	$(FC) -o $@  $(source)

.PHONY : clean

clean:
	@rm -f *.o *.mod rand
