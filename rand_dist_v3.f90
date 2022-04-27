!include "nrtype.f90"
!include "nrutil.f90"
!include "spline_mod.f90"

! This module is useful for generating random numbers which follows any kind of arbitrary 
! probability distribution function (pdf)
! Initiate the random generator by first calling the init_generator subroutine, the given
! pdf should be named as 'xs_dist.txt', gen_randist subroutine then generated random numbers
! the main program should always end with stop_generator subroutine
! Author: Abhijit Baishya
! Date: 27.04.2022

module rand_dist
use spline_mod
implicit none

private
real,allocatable :: x(:),p(:),y2(:),ji(:),sol2(:),sol(:)
real :: u,res
integer :: N 

public :: init_generator,stop_generator,gen_randist

contains



! Initiating the random generator
subroutine init_generator(fname) 
implicit none 

character(*), intent(in) :: fname
integer :: i,j,iostat

open(10,file=fname) 

N = 0
do 
    read(10,*,iostat=iostat)
    if (iostat < 0) exit
    N = N+1
enddo

rewind(10) 

allocate(x(N),p(N),y2(N),ji(N),sol(N),sol2(N))

do i=1,N
    read(10,*)x(i),p(i) 
enddo 

call spline(x,p,0.,0.,y2) 
call simpsint(x(1),x(N),N,f,res)
!print*,res
p = p/res 
call spline(x,p,0.,0.,y2)

do i=1,N
    ji(i) = CDF(x(i))
    call solve_cdf(ji(i),sol(i))
    !print*,ji(i),sol(i)
enddo 

call spline(ji,sol,0.,0.,sol2) 

close(10) 

return
end subroutine


! Stopping the random generator
subroutine stop_generator()
implicit none 
deallocate(x,p,y2,ji,sol,sol2) 
return 
end subroutine 



! Interpolated function
function f(xx)
use spline_mod
implicit none 
real,intent(in) :: xx
real :: f
f = splint(x,p,y2,xx)
end function 



!!Simpson's 1/3 rule integration
subroutine simpsint(a, b, N, f, integral)
implicit none 
real, external :: f 
real,intent(in) :: a,b
real,intent(out) :: integral
integer,intent(in) :: N
real :: h 
integer :: i

h = (b-a)/N
integral = 0.0

do i=1,N
integral = integral + (h/6.0)*( f(a+(i*1.0-1.0)*h) + 4.0*f(a+(i*1.0-.5)*h) + f(a+i*h) )
enddo

return
end subroutine


!Cumulative distribution function
function CDF(xx)
implicit none
!real, external :: f
!real,intent(in) :: x(:),p(:),xx 
real :: a,b,step,CDF,xx
integer :: N 
a = minval(x) ; b = xx 
step = x(2)-x(1) ! For equal steps 
N = nint((b-a)/step)
call simpsint(a,b,N,f,CDF)
return 
end function


! Solving Inverse CDF to find random number following distribution
subroutine solve_cdf(y,sol) 
use spline_mod
implicit none 
!real, external :: CDF 
real,intent(in) :: y 
real,intent(out) :: sol
real :: min,max,step 
real,allocatable :: soln(:),testx(:),testp(:) 
integer :: i,j,N,locmin(1),locmax(1)
N = size(x) 
!step = x(2)-x(1) 
allocate(soln(N),testx(N),testp(N))
testx = x ; testp = p
min = minval(x) ; max = maxval(x)
do j=1,1000
    do i=1,N 
        soln(i) = abs(y-CDF(testx(i)))
    enddo  
    
    locmin = minloc(soln) ; locmax = maxloc(soln)
    sol = testx(locmin(1))
    !print*,sol,y,CDF(sol)
    if (abs(y-CDF(sol)) < 1.0e-06) exit 
    min = testx(locmin(1)-1) ; max = testx(locmin(1)+1)
    step = (max-min)/N
    do i=1,N 
        testx(i) = min + (i-1)*step
        testp(i) = splint(x,p,y2,testx(i))
    enddo
enddo 
deallocate(soln,testx,testp)
return 
end subroutine 


! Random Number Generator following distribution
subroutine gen_randist(r) 
implicit none 
real,intent(out) :: r
integer :: j
call random_number(u) 
do j=1,N 
    if (u > ji(j) .and. u < ji(j+1)) then 
        r = splint(ji,sol,sol2,u)
    endif
enddo
return 
end subroutine 

end module
