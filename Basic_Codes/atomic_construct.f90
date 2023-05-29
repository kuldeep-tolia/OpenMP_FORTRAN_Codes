! Works as an alternative for critical block when "atomic" operations are used

integer function increase(thread_id, j)		!Reports thread id performing addition operation to s2 variable

implicit none

    integer, intent(in) :: thread_id, j
    
    print *, "Function run by thread number =", thread_id
    increase = j

end function increase

program atomic_construct

use omp_lib
implicit none
    
    integer, parameter :: num_threads = 4, m = 20
    integer :: thread_id, i, j, s1 = 0, s2 = 0
    integer, external :: increase
    
    thread_id = 0
    !$ call omp_set_num_threads(num_threads)
    
    !$omp parallel do private(thread_id, j) shared(s1, s2)
    	
    	do i = 0, m-1
            j = 2 * i - i
            !$ thread_id = omp_get_thread_num()
            !$omp atomic
                s1 = s1 + j
                s2 = s2 + increase(thread_id, j)
        end do
    
    !$omp end parallel do
    
    print *, "checking sum-1 = ", s1
    print *, "checking sum-2 = ", s2

end program atomic_construct
