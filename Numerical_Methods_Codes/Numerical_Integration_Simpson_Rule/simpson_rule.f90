program tut2q3

use omp_lib
implicit none

	double precision :: a, b, integration_result, error
	double precision, parameter :: exact_result = 0.1985730d0
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
	
	n = 32				! number of divisions for integration
	a = 1.0d0			! integration lower limit
	b = 4.0d0 * atan(1.0d0)		! integration upper limit
	integration_result = 0.0d0
		
	call simp_rule(a, b, n, integration_result)
	
	error = abs(exact_result - integration_result)
	write(*, 100) "Number of threads =", num_thread, " Number of Simpson divisions =", n
	100 format (A, I3, A, I5)
	write(*, 101) "The integration for the given function between limits 1.0 and PI =", integration_result
	101 format (A, F15.7)
	write(*, 102) "Absolute error with the exact integration value =", error
	102 format (A, F15.7)
	
contains

	function func(x) result(y)
	
	implicit none
		
		double precision :: x, y
		
		y = sin(x) / (2.0d0 * x**3)
		
	end function func
	
	subroutine simp_rule(x_s, x_e, n, global_result)
	
	implicit none
		
		double precision :: x_s, x_e, global_result
		integer :: n
		double precision :: h, x, local_sum
		integer :: i
		
		h = (x_e - x_s) / n
		local_sum = func(x_s) + func(x_e)
		
		!$omp parallel num_threads(num_thread) reduction(+:local_sum)
			
			!$omp do private(i)
			do i = 1, n-1, 2
				local_sum = local_sum + 4.0d0 * func(x_s + i * h)
			end do
			!$omp end do
			
			!$omp do private(i)
			do i = 2, n-2, 2
				local_sum = local_sum + 2.0d0 * func(x_s + i * h)
			end do
			!$omp end do	
					
		!$omp end parallel
		
		global_result = local_sum * h / 3.0d0
		
	end subroutine simp_rule
	
end program tut2q3
