# randDist
Written using object oriented FORTRAN

Generate random numbers following uniform, normal as well as any arbitrary probability distribution function
class methods:

i) subroutine uni(randomNumber, lowerLimit, upperLimit): Generate uniformly distributed random number (integer, real or double precision) or randomNumber array, lower and upper limits are optional and if not given will generate numbers between 0 and 1

ii) subroutine normal(randomNumber, mean, sigma): Generate normally distributed random number (integer, real or double precision) or randomNumber array

iii) subroutine sphere(theta, phi): Generate uniformly distributed polar angles (array also possible) inside a sphere 

iv) subroutine rndm(randomNumber): Generate arbitrarily distributed (distribution has to be provided via a file) random number (integer, real or double precision) or randomNumber array; two helper subroutines are required for this, 1. subroutine init(filename) and subroutine finish() 

gfortran compiler is needed.

test_rand.f90 is a sample main program. 

To run the smaple main program run "make" in the terminal; an executable "rand" will be created, otherwise barring the test_rand.f90 file other files are necessary as the module rand needs nrtype.f90, nrutil.f90 and spline_mod.f90 .

The sample probability distribution function is given in xs_dist.txt file.
