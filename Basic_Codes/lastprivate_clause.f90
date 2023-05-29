! lastprivate clause:
    !   Each thread has its own copy of variables to modify
    !   Values are deleted after the scope of clause
    !   private variables will be auto-initialized to 0 before usage!
    !   after scope, the variables assume the last value inside the scope
    
program last_private_clause

use omp_lib
implicit none

	integer, parameter :: num_threads = 4
	integer :: i, thread_id
	real*8 :: t1, t2, ep
	
	!$ call omp_set_num_threads(num_threads)
	call cpu_time(t1)
	
	print*, "Thread value before lastprivate: ", thread_id
	print*, "Parallel Mode ON!!!"
	!$omp parallel do lastprivate(thread_id)
		do i = 1, 12
			thread_id = omp_get_thread_num() + 20 + i
			print*, "i = ", i, " thread = ", omp_get_thread_num(), " thread value = ", thread_id
		end do
	!$omp end parallel do
	call cpu_time(t2)
	ep = (t2 - t1) / num_threads
	print*, "Thread value after lastprivate = ", thread_id
	print*, "Computing time = ", ep, " s"

end program last_private_clause
