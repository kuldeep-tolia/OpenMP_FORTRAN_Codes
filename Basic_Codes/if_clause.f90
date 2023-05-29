!Program to demonstrate if-clause in OpenMP

program if_clause

use omp_lib
implicit none

	integer, parameter :: thread_num = 2
	integer :: thread_id = -1
	
	!$ call omp_set_num_threads(thread_num)
	!!$omp parallel if (thread_num > 1) private(thread_id)
       !$omp parallel if (thread_num > 1) firstprivate(thread_id)
       	
       	!$thread_id = omp_get_thread_num()
       	print*, "if clause executed by thread id =", thread_id
       	
       !$omp end parallel
       
       print*, "value of thread id after parallel block =", thread_id

end program if_clause
