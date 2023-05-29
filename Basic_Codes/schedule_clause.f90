program schedule_clause

use omp_lib
implicit none
    
    integer, parameter :: num_threads = 4, m = 15
    integer :: thread_id, i
    
    thread_id = 0
    !$ call omp_set_num_threads(num_threads)
    
    !$omp parallel do private(thread_id) schedule(guided, 4)
        do i = 0, m-1
            !$ thread_id = omp_get_thread_num()
            print *, "Iteration ", i, " done by thread ", thread_id
        end do
    !$omp end parallel do

end program schedule_clause
