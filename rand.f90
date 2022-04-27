! SAMPLE TEST MAIN PROGRAM
program rand_test 
use rand_dist 
implicit none 

real :: r
integer :: i,j,N_event

N_event = 100000
open(11,file="xs_algo.txt") 

call init_generator("xs_dist.txt") 

write(*,*)"--------------------100%"
do i=1,N_event
    if (modulo(i*(100.d0/N_event),5.d0) == 0.d0) then 
        write(*,'(a1)',advance='no')"#"
    endif
    call gen_randist(r)
    write(11,*)r
enddo 
write(*,'(a6)')'#DONE!'
call stop_generator() 
close(11)
end program