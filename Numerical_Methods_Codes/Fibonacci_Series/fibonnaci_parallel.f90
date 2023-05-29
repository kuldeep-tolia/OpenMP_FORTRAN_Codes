program fibonacci_parallel

use omp_lib
implicit none

	integer, parameter :: N = 20
	integer :: k, my_id
	integer :: num_thread = 1
	character(100) :: numchar
	character(100) :: nameq
	
	if(COMMAND_ARGUMENT_COUNT() .ne. 1) then
		write(*, *) "Command line argument of number of threads is required."
		STOP
	end if
	
	call GET_COMMAND_ARGUMENT(0, nameq)	
	call GET_COMMAND_ARGUMENT(1, numchar)
	read(numchar, *) num_thread
	
	!$omp parallel num_threads(num_thread) private(my_id)
		!$omp single
			my_id = omp_get_thread_num()
			write(*, *) "N     fib(N)"
			do k = 0, N
				write(*, *) k, fib(k)
			end do
		!$omp end single
	!$omp end parallel

contains
	
	recursive function fib(n) result(res)
	
	implicit none
		
		integer :: n, res, i, j
			
		if(n .lt. 2) then
			res = n
		else
			!$omp task shared(i)
				i = fib(n-1)
			!$omp end task
			
			!$omp task shared(j)
				j = fib(n-2)
			!$omp end task
			
			!$omp taskwait
			res = i + j
		end if
	
	end function fib	

end program fibonacci_parallel
