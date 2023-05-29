program nested_parallel

use omp_lib
implicit none

	integer:: id
	call omp_set_nested(.true.)
	call omp_set_num_threads(2)
	!$omp parallel private(id)
		id = omp_get_thread_num()
		write(*, *) "Hello from thread=", id
	!$omp parallel private(id)
		id = omp_get_thread_num()
		write(*, *) "Another HI from thread=", id
	!$omp end parallel
	!$omp end parallel

end program nested_parallel
