include "nrtype.f90"
include "nrutil.f90"
include "spline_mod.f90"

program rand_dist 
use spline_mod
implicit none

real,allocatable :: x(:),p(:),y2(:)
real :: u,res,z,r 
integer :: i,j,iostat,N,N_event 


open(10,file="xs_dist_exp.txt") 
open(11,file="xs_algo_exp.txt") 

N = 0
do 
    read(10,*,iostat=iostat)
    if (iostat < 0) exit
    N = N+1
enddo

rewind(10) 

allocate(x(N),p(N),y2(N))

do i=1,N
    !call random_number(u) 
    !x = sqrt(-2.0*log(u))
    !x(i) = -2.0 + i*(8.0/N)
    !write(10,*)x(i),exp(-(x(i)-2.0)**2.0)
    read(10,*)x(i),p(i) 
enddo 
!print*,p(1:10)
call spline(x,p,0.,0.,y2) 
!print*,y2(1:5)!splint(x,p,y2,60.5)


call simpsint(x(1),x(N),N,f,res)
!print*,res
p = p/res 
call spline(x,p,0.,0.,y2)
!call simpsint(x(1),x(N),N,f,res)
!print*,res!,CDF(x,p,f,5.0)
!do i=1,5
!call random_number(u) 
!print*,u 
!enddo
!print*,"uniform random number",u
!call solve_cdf(u,r)
!print*,"dist. random number",r

N_event = 100000
j = 0
write(*,*)"--------------------100%"
do i=1,N_event
    if (modulo(i*(100.d0/N_event),5.d0) == 0.d0) then 
        j = j+1
        !print*,"completed",j*5.0,"% at event no:", i
        write(*,'(a1)',advance='no')"#"
    endif
    call random_number(u) 
    !x = sqrt(-2.0*log(u))
    !z = -2.0 + i*(8.0/(2*N))
    !write(10,*)x(i),exp(-(x(i)-2.0)**2.0)
    call solve_cdf(u,r)
    write(11,*)r
    !write(1,*)z,splint(x,p,y2,z)
enddo 
write(*,'(a6)')'#DONE!'

deallocate(x,p,y2) 

close(10) 
close(11)

contains 


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



end program
