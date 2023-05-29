! OpenMP program to demonstrate the critical clause    
program critical_clause

use omp_lib
implicit none

	integer, parameter :: num_threads = 4
	integer :: i, thread_id
	real*8 :: t1, t2
	
	!$ call omp_set_num_threads(num_threads)
	t1 = omp_get_wtime()
	
	print*, "Thread value before lastprivate: ", thread_id
	print*, "Parallel Mode ON!!!"
	!$omp parallel do lastprivate(thread_id)
		do i = 1, 12
			!$omp critical
			thread_id = omp_get_thread_num() + 20 + i
			print*, "i = ", i, " thread = ", omp_get_thread_num(), " thread value = ", thread_id
			!$omp end critical
		end do
	!$omp end parallel do
	t2 = omp_get_wtime()
	print*, "Thread value after lastprivate = ", thread_id
	print*, "Computing time = ", t2-t1, " s"

end program critical_clause
