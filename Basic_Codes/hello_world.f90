program Hello_World

use omp_lib
implicit none

	integer :: id, num_threads = 1
	character(100) :: numchar
	character(100) :: nameq
	
	if(COMMAND_ARGUMENT_COUNT() .ne. 1) then
		write(*, *) "Command line argument of number of threads is required."
		STOP
	end if
	
	call GET_COMMAND_ARGUMENT(0, nameq)
	write(*, *) nameq
	
	call GET_COMMAND_ARGUMENT(1, numchar)
	read(numchar, *) num_threads
	
	!$ call omp_set_num_threads(num_threads)
	
	!$omp parallel private(id)
!	!$omp critical				!identifies a section of code that must be executed by a single thread at a time. 

		id = omp_get_thread_num()
		write(*, *) "Hello World from thread id =", id, " out of total threads = ", num_threads
		
!	!$omp end critical	
	!$omp end parallel 

end program Hello_World
