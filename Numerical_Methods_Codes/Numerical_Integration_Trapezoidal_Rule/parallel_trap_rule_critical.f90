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
	
	!$omp parallel num_threads(num_thread) 
		
		call trap_rule(a, b, n, integration_result)
		
	!$omp end parallel
	
	write(*, 100) "The integration for the given function between limits", a, " and", b, " =", integration_result
	100 format (A, F10.2, A, F10.2, A, F15.7)
	
contains
	
	subroutine trap_rule(x_s, x_e, n, global_result)
	
	implicit none
		
		double precision :: x_s, x_e, global_result
		integer :: n
		double precision :: h, x, local_sum, local_a, local_b
		integer :: i, my_id, local_n
		
		h = (x_e - x_s) / n
		local_n = n / num_thread
		local_a = x_s
		local_b = x_e
		my_id = omp_get_thread_num()
		local_a = x_s + my_id * local_n * h
		local_b = local_a + local_n * h
		local_sum = ((1.0d0 + sin(local_a)) + (1.0d0 + sin(local_b))) / 2.0d0	! given function to integrate is 1.0 + sin(x)
		
		do i = 1, local_n-1
			x = local_a + i * h
			local_sum = local_sum + 1.0d0 + sin(x)
		end do
		
		local_sum = local_sum * h
		
		! to avoid race condition
		!$omp critical		
			global_result = global_result + local_sum
		!$omp end critical
		
	end subroutine trap_rule
	
end program parallel_trapezoidal_rule
