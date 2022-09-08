! PROGRAM TO TEST THE RANDOM MODULE
program test_random 
use rand
implicit none 

type (random) :: r 
real :: a,a1(5),a2,a3(5) 
real*8 :: b,b1(5),b2,b3(5) 
integer :: i,j,N_event,k,k1(5) 

! CALLING USUAL UNIFORM RANDOM NUMBERS
call r%uniform(a)
call r%uniform(b)
!print*,a,b
! CALLING UNIFORM RANDOM NUMBERS WITH LOWER AND UPPER LIMIT
call r%uniform(a,2.0,10.0)
call r%uniform(b,10.d0,15.d0)
!print*,a,b
! CALLING USUAL UNIFORM INTGER RANDOM NUMBERS
call r%uniform(k,5,10)
!print*,k
! CALLING USUAL UNIFORM RANDOM NUMBER ARRAY
call r%uniform(a1,2.0,15.0)
call r%uniform(b1,2.d0,15.d0)
!print*,a1,b1
! CALLING USUAL UNIFORM INTEGER RANDOM NUMBER ARRAY
call r%uniform(k1,5,10)
!print*,k1
! CALLING NORMAL RANDOM NUMBERS
call r%normal(a,5.0,2.0)
call r%normal(b,5.d0,2.d0)
!print*,a,b
! CALLING NORMAL RANDOM NUMBERS ARRAY
call r%normal(a1,5.0,2.0)
call r%normal(b1,5.d0,2.d0)
!print*,a1,b1

N_event = 100000
open(11,file="xs_algo.txt") 

call r%init("xs_dist.txt") 

write(*,*)"--------------------100%"
do i=1,N_event
    if (modulo(i*(100.d0/N_event),5.d0) == 0.d0) then 
        write(*,'(a1)',advance='no')"#"
    endif
    call r%rndm(b)
    write(11,*)b
enddo 
write(*,'(a6)')'#DONE!'
call r%rndm(a1)
call r%rndm(b1)
print*,a1,b1
call r%finish() 
close(11)

call r%sphere(a,a2)
call r%sphere(b,b2) 
call r%sphere(a1,a3)
call r%sphere(b1,b3) 
!print*,a,a2,b,b2,a1,a3,b1,b3  

end program