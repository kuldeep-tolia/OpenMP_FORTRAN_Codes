! OpenMP parallelized program to calculate the value of PI using series expansion method
program calc_pi

use omp_lib
implicit none

	integer, parameter :: N = 100000000
	integer :: num_thread, k
	double precision :: sig = 1.0d0, summation = 0.0d0, pi = 0.0d0, t_start, t_end
	character(100) :: numchar
	character(100) :: nameq
	
	if(COMMAND_ARGUMENT_COUNT() .ne. 1) then
		write(*, *) "Command line argument of number of threads is required."
		STOP
	end if
	
	call GET_COMMAND_ARGUMENT(0, nameq)	
	call GET_COMMAND_ARGUMENT(1, numchar)
	read(numchar, *) num_thread
	
	t_start = omp_get_wtime()
	!$omp parallel num_threads(num_thread) reduction(+:summation) default(shared) private(k, sig)
		!$omp do
		do k = 0, N-1
			if(mod(k, 2) .eq. 0) then
				sig = 1.0d0
			else
				sig = -1.0d0
			end if
			
			summation = summation + (sig / (2.0d0 * k + 1.0d0))
		end do
		!$omp end do
	!$omp end parallel
	t_end = omp_get_wtime()
	
	pi = 4.0d0 * summation
	write(*, *) "Approximate value of PI using series expansion is ", pi
	write(*, *) "Execution time = ", t_end-t_start, "seconds"	

end program calc_pi
