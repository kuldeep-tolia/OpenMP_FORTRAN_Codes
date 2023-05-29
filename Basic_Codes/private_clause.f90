! Private clause:
    !   Each thread has its own copy of variables to modify
    !   Values are deleted after the scope of clause
    !   private variables are auto-initialized to 0 before usage!
    !   after scope, the variables assume the value prior to the scope

program private_clause

use omp_lib
implicit none

	integer, parameter :: num_threads = 4
	integer :: i, thread_id = 11
	
	!$ call omp_set_num_threads(num_threads)
	
	print*, "Thread value before private: ", thread_id
	print*, "Parallel Mode ON!!!"
	!$omp parallel private(thread_id)
		print*, "Thread value = ", thread_id
		thread_id = omp_get_thread_num()
		print*, "Active thread = ", thread_id
		thread_id = thread_id + 12
		print*, "Thread = ", omp_get_thread_num(), " and thread value = ", thread_id
	!$omp end parallel
	print*, "Thread value after private = ", thread_id

end program private_clause
