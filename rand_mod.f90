! Random number genrator module
! Author: Abhijit Baishya 
! Last Edited: 08/09/2022

module rand
use spline_mod
implicit none 

type, public :: random
    private
    real*8,allocatable :: x(:),p(:),y2(:),ji(:),sol2(:),sol(:)
    real*8 :: res
    integer :: N 

    contains 

        ! Uniform distributions
        procedure, pass, private :: uni
        procedure, pass, private :: uni_d
        procedure, pass, private :: uni_lim
        procedure, pass, private :: uni_lim_d
        procedure, pass, private :: uni_i
        procedure, pass, private :: uni_arr
        procedure, pass, private :: uni_arr_d
        procedure, pass, private :: uni_lim_arr
        procedure, pass, private :: uni_lim_arr_d
        procedure, pass, private :: uni_arr_i

        ! NORMAL DISTRIBUTION 
        procedure, pass, private :: norm 
        procedure, pass, private :: norm_d 
        procedure, pass, private :: norm_arr 
        procedure, pass, private :: norm_arr_d 

        ! UNIFORMLY DISTRIBUTED RANDOM SPHERICAL ANGLES FOR A SPHERE 
        procedure, pass, private :: rand_sphr_sp
        procedure, pass, private :: rand_sphr_dp 
        procedure, pass, private :: rand_sphr_arr_sp
        procedure, pass, private :: rand_sphr_arr_dp

        ! ARBITRARY DISTRIBUTION 
        procedure, pass, private :: init_generator
        procedure, pass, private :: stop_generator
        procedure, pass, private :: gen_randist 
        procedure, pass, private :: gen_randist_arr
        procedure, pass, private :: gen_randist_d
        procedure, pass, private :: gen_randist_arr_d
        procedure, pass, private :: f 
        procedure, pass, private :: simpsint
        procedure, pass, private :: CDF 
        procedure, pass, private :: solve_cdf
        
        

        ! PUBLIC PROCEDURES (AKA CLASS METHODS)
        ! 'uniform' subroutine will generate uniformly distributed random numbers 
        generic, public :: uniform => uni,uni_d,uni_lim,uni_lim_d,uni_i,uni_arr,uni_arr_d,&
        uni_lim_arr,uni_lim_arr_d,uni_arr_i

        ! 'normal' subroutine will generate uniformly distributed random numbers 
        generic, public :: normal => norm,norm_d,norm_arr,norm_arr_d

        ! 'sphere' subroutine will generate uniformly distributed random polar angles 
        generic, public :: sphere => rand_sphr_sp,rand_sphr_dp,rand_sphr_arr_sp,rand_sphr_arr_dp

        ! 'init', 'finish' & 'rndm' subroutine will generate random numbers following a arbitrary 
        ! distribution from a file name given in the 'init' subroutine
        generic, public :: init => init_generator
        generic, public :: finish => stop_generator
        generic, public :: rndm => gen_randist,gen_randist_d,gen_randist_arr,gen_randist_arr_d

end type 

contains 

! UNIFORM DISTRIBUTIONS 
subroutine uni(this,r) ! Single precision
    class(random) :: this 
    real :: r 
    call random_number(r)
    return 
end subroutine 

subroutine uni_d(this,r)  ! Double precision
    class(random) :: this 
    real*8 :: r 
    call random_number(r)
    return 
end subroutine 

subroutine uni_lim(this,r,low,high) ! Single precision between limits
    class(random) :: this 
    real :: r,low,high 

    call random_number(r)
    r = low + (high-low)*r 
    return 
end subroutine 

subroutine uni_lim_d(this,r,low,high) ! Double precision between limits
    class(random) :: this 
    real*8 :: r,low,high 
    
    call random_number(r)
    r = low + (high-low)*r 
    return 
end subroutine 

subroutine uni_i(this,r,low,high) ! Random integers between limits
    class(random) :: this 
    integer :: r,low,high 
    real :: s 
    call random_number(s)
    r = low + int((high-low)*s) 
    return 
end subroutine 

subroutine uni_arr(this,r) ! single precision random number array of arbitrary size
    class(random) :: this 
    real :: r(:) 
    call random_number(r)
    return 
end subroutine 

subroutine uni_arr_d(this,r) ! double precision random number array of arbitrary size
    class(random) :: this 
    real*8 :: r(:) 
    call random_number(r)
    return 
end subroutine 

! single precision random number array of arbitrary size between limits
subroutine uni_lim_arr(this,r,low,high) 
    class(random) :: this 
    real :: r(:),low,high 

    call random_number(r)
    r = low + (high-low)*r 
    return 
end subroutine 

! double precision random number array of arbitrary size between limits
subroutine uni_lim_arr_d(this,r,low,high)
    class(random) :: this 
    real*8 :: r(:),low,high 

    call random_number(r)
    r = low + (high-low)*r 
    return 
end subroutine 

! integer random number array of arbitrary size between limits
subroutine uni_arr_i(this,r,low,high)
    class(random) :: this 
    integer :: r(:),low,high 
    real,allocatable :: s(:)
    integer :: n 
    
    n = size(r)
    allocate(s(n))
    call random_number(s)
    r = low + int((high-low)*s) 
    deallocate(s)
    return 
end subroutine 
! END OF UNIFORM DISTRIBUTIONS

! NORMAL DISTRIBUTIONS 
! Normally distributed single precision random numbers with mean and sigma 
subroutine norm(this,r,mean,sigma)
    class(random) :: this 
    real :: mean,sigma
    real :: r,r1,r2
    real, parameter :: pi=4.0*atan(1.0)
    call random_number(r1)
    call random_number(r2)
    r = sqrt(-2.0*log(r1))*cos(2.0*pi*r2)
    r = mean + sigma*r
    return
end subroutine norm

! Normally distributed double precision random numbers with mean and sigma
subroutine norm_d(this,r,mean,sigma)
    class(random) :: this 
    real*8 :: mean,sigma
    real*8 :: r,r1,r2
    real*8, parameter :: pi=4.d0*atan(1.d0)
    call random_number(r1)
    call random_number(r2)
    r = sqrt(-2.0*log(r1))*cos(2.0*pi*r2)
    r = mean + sigma*r
    return
end subroutine norm_d

! Normally distributed single precision random number array with mean and sigma
subroutine norm_arr(this,r,mean,sigma)
    class(random) :: this 
    real :: mean,sigma
    real :: r(:)
    real,allocatable :: r1(:),r2(:)
    real, parameter :: pi=4.0*atan(1.0)
    integer :: n 
    n = size(r) 
    allocate(r1(n),r2(n)) 
    call random_number(r1)
    call random_number(r2)
    r = sqrt(-2.0*log(r1))*cos(2.0*pi*r2)
    r = mean + sigma*r
    deallocate(r1,r2)
    return
end subroutine norm_arr

! Normally distributed double precision random number array with mean and sigma
subroutine norm_arr_d(this,r,mean,sigma)
    class(random) :: this 
    real*8 :: mean,sigma
    real*8 :: r(:)
    real*8,allocatable :: r1(:),r2(:)
    real*8, parameter :: pi=4.d0*atan(1.d0)
    integer :: n 
    n = size(r) 
    allocate(r1(n),r2(n)) 
    call random_number(r1)
    call random_number(r2)
    r = sqrt(-2.0*log(r1))*cos(2.0*pi*r2)
    r = mean + sigma*r
    deallocate(r1,r2)
    return
end subroutine norm_arr_d
! END OF NORMAL DISTRIBUTIONS 


! UNIFORMLY GENERATED VECTORS INSIDE A SPHERE
! GENERATES THETA,PHI SUCH THAT UNIFORM DISTRIBUTION ON SPHERE IS POSSIBLE
subroutine rand_sphr_sp(this,theta,phi) 
class(random) :: this
real,intent(out) :: theta,phi
real :: r,s 
real, parameter :: pi=4.0*atan(1.0), d2r=pi/180.0, r2d=1.0/d2r
call random_number(r)
call random_number(s) 
theta = (acos(1.0-2.0*r))*r2d
phi = (2.0*pi*s)*r2d
return 
end subroutine

! DOUBLE PRECISION VERSION OF THE PREVIOUS SUBROUTINE
subroutine rand_sphr_dp(this,theta,phi) 
class(random) :: this
real*8,intent(out) :: theta,phi
real*8 :: r,s 
real*8, parameter :: pi=4.0*atan(1.d0), d2r=pi/180.d0, r2d=1.d0/d2r
call random_number(r)
call random_number(s) 
theta = (acos(1.0-2.0*r))*r2d
phi = (2.0*pi*s)*r2d
return 
end subroutine

! GENERATES THETA,PHI ARRAY OF ARBITRARY SIZE SUCH THAT UNIFORM DISTRIBUTION ON SPHERE IS POSSIBLE
subroutine rand_sphr_arr_sp(this,theta,phi) 
class(random) :: this
real,intent(out) :: theta(:),phi(:)
real,allocatable :: r(:),s(:) 
real, parameter :: pi=4.d0*atan(1.0), d2r=pi/180.0, r2d=1.0/d2r
integer :: n 
n = size(theta)
allocate(r(n),s(n)) 
call random_number(r)
call random_number(s) 
theta = (acos(1.0-2.0*r))*r2d
phi = (2.0*pi*s)*r2d
deallocate(r,s)
return 
end subroutine

! DOUBLE PRECISION VERSION OF THE PREVIOUS SUBROUTINE
subroutine rand_sphr_arr_dp(this,theta,phi) 
class(random) :: this
real*8,intent(out) :: theta(:),phi(:)
real,allocatable :: r(:),s(:) 
real*8, parameter :: pi=4.d0*atan(1.d0), d2r=pi/180.d0, r2d=1.d0/d2r
integer :: n 
n = size(theta)
allocate(r(n),s(n)) 
call random_number(r)
call random_number(s) 
theta = (acos(1.0-2.0*r))*r2d
phi = (2.0*pi*s)*r2d
deallocate(r,s)
return 
end subroutine

! CUSTOM DISTRIBUTION
! InitiaL routine for the arbitrary random generator, file name of the distribution 
! function in (x,y) format has be given as an input argument
subroutine init_generator(this,fname)  
    class(random) :: this 
    character(*), intent(in) :: fname
    integer :: i,j,iostat

    open(10,file=fname) 

    ! DETERMING THE DATA SIZE
    this%N = 0
    do 
        read(10,*,iostat=iostat)
        if (iostat < 0) exit
        this%N = this%N+1
    enddo
    
    rewind(10) 

    if (allocated(this%x)) deallocate(this%x,this%p,this%y2,this%ji,this%sol,this%sol2)
    allocate(this%x(this%N),this%p(this%N),this%y2(this%N),this%ji(this%N),this%sol(this%N),this%sol2(this%N))

    ! READING DISTRIBUTION FROM FILE
    do i=1,this%N
        read(10,*)this%x(i),this%p(i) 
    enddo 

    call spline(this%x,this%p,0.d0,0.d0,this%y2) 
    call this%simpsint(this%x(1),this%x(this%N),this%N,this%res)
    this%p = this%p/this%res 

    call spline(this%x,this%p,0.d0,0.d0,this%y2)

    ! SOLVING CUMULATIVE DISTRIBUTION FUNCTION
    do i=1,this%N
        this%ji(i) = this%CDF(this%x(i))
        call this%solve_cdf(this%ji(i),this%sol(i))
    enddo 

    ! SAVING INTERPOLATED SOLUTIONS IN SOL ARRAY
    call spline(this%ji,this%sol,0.d0,0.d0,this%sol2) 

    close(10) 
    return
end subroutine


! Stopping the random generator
subroutine stop_generator(this)
    class(random) :: this
    deallocate(this%x,this%p,this%y2,this%ji,this%sol,this%sol2) 
    return 
end subroutine 

! Interpolated function
function f(this,xx)
    class(random) :: this
    real*8,intent(in) :: xx
    real*8 :: f
    f = splint(this%x,this%p,this%y2,xx)
    return 
end function 

!!Simpson's 1/3 rule integration
subroutine simpsint(this,a, b, N, integral)
    class(random) :: this
    real*8,intent(in) :: a,b
    real*8,intent(out) :: integral
    integer,intent(in) :: N
    real*8 :: h 
    integer :: i

    h = (b-a)/N
    integral = 0.0

    do i=1,N
    integral = integral + (h/6.0)*( this%f(a+(i*1.0-1.0)*h) + 4.0*this%f(a+(i*1.0-.5)*h) + this%f(a+i*h) )
    enddo

    return
end subroutine


!Cumulative distribution function
function CDF(this,xx)
    class(random) :: this
    real*8 :: a,b,step,CDF,xx
    integer :: N 
    a = minval(this%x) ; b = xx 
    step = this%x(2)-this%x(1) ! For equal steps 
    N = nint((b-a)/step)
    call this%simpsint(a,b,N,CDF)
    !print*,a,b,N,CDF
    return 
end function


! Solving Inverse CDF to find random number following distribution
subroutine solve_cdf(this,y,sol) 
    class(random) :: this
    real*8,intent(in) :: y 
    real*8,intent(out) :: sol
    real*8 :: min,max,step 
    real*8,allocatable :: soln(:),testx(:),testp(:) 
    integer :: i,j,N,locmin(1),locmax(1)
    N = size(this%x) 
    allocate(soln(N),testx(N),testp(N))
    testx = this%x ; testp = this%p
    min = minval(this%x) ; max = maxval(this%x)
    do j=1,1000
        do i=1,N 
            soln(i) = abs(y-this%CDF(testx(i)))
        enddo  

        locmin = minloc(soln) ; locmax = maxloc(soln)
        sol = testx(locmin(1))
        if (abs(y-this%CDF(sol)) < 1.0e-06) exit 
        min = testx(locmin(1)-1) ; max = testx(locmin(1)+1)
        step = (max-min)/N
        do i=1,N 
            testx(i) = min + (i-1)*step
            testp(i) = splint(this%x,this%p,this%y2,testx(i))
        enddo
    enddo 
    deallocate(soln,testx,testp)
    return 
end subroutine 


! generates random number following distribution, single precision version
subroutine gen_randist(this,r) 
    class(random) :: this
    real,intent(out) :: r
    real*8 :: u 
    call random_number(u) 
    r = splint(this%ji,this%sol,this%sol2,u)
    return 
end subroutine 


! generates random number following distribution, double precision version
subroutine gen_randist_d(this,r) 
    class(random) :: this
    real*8,intent(out) :: r
    real*8 :: u 
    call random_number(u)  
    r = splint(this%ji,this%sol,this%sol2,u)
    return 
end subroutine


! generates random number array following distribution, single precision version
subroutine gen_randist_arr(this,r) 
    class(random) :: this
    real,intent(out) :: r(:)
    real*8 :: u 
    integer :: i,n 
    n = size(r) 
    do i=1,n 
        call random_number(u) 
        r(i) = splint(this%ji,this%sol,this%sol2,u)
    enddo
    return 
end subroutine 


! generates random number array following distribution, double precision version
subroutine gen_randist_arr_d(this,r) 
    class(random) :: this
    real*8,intent(out) :: r(:)
    real*8 :: u 
    integer :: i,n 
    n = size(r) 
    do i=1,n 
        call random_number(u) 
        r(i) = splint(this%ji,this%sol,this%sol2,u)
    enddo
    return 
end subroutine


end module rand


