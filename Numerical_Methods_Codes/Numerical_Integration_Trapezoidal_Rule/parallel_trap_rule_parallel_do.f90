program parallel_trapezoidal_rule

use omp_lib
implicit none

	double precision :: a, b, integration_result
	integer :: n
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
	
	n = 72				! number of divisions for integration
	a = 0				! integration lower limit
	b = 4.0d0 * atan(1.0d0)		! integration upper limit
	integration_result = 0.0d0
		
	call trap_rule(a, b, n, integration_result)
		
	write(*, 100) "The integration for the given function between limits", a, " and", b, " =", integration_result
	100 format (A, F10.2, A, F10.2, A, F15.7)
	
contains
	
	subroutine trap_rule(x_s, x_e, n, global_result)
	
	implicit none
		
		double precision :: x_s, x_e, global_result
		integer :: n
		double precision :: h, x
		integer :: i, thread_count
		
		thread_count = omp_get_num_threads()		
		h = (x_e - x_s) / n
		global_result = ((1.0d0 + sin(x_s)) + (1.0d0 + sin(x_e))) / 2.0d0	! given function to integrate is 1.0 + sin(x)
		
		!$omp parallel num_threads(thread_count) reduction(+:global_result)
			!$omp do private(i)
				
				do i = 1, n-1
					global_result = global_result + 1.0d0 + sin(x_s + i * h)
				end do
				
			!$omp end do
		!$omp end parallel
			
		global_result = global_result * h
		
	end subroutine trap_rule
	
end program parallel_trapezoidal_rule
