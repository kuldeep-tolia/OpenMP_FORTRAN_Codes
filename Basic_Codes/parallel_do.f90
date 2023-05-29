program parallel_do_execution

use omp_lib
implicit none

	integer :: num_threads = 2
	integer, parameter :: m = 100, n = 100000, o = 100
	integer :: i, j, k
	real*8, dimension(1:m, 1:n) :: a
	real*8, dimension(1:n, 1:o) :: b
	real*8, dimension(1:m, 1:o) :: c
	real*8 :: t1, t2, ep
	
	a(:, :) = 1.0
	b(:, :) = 1.0
	c(:, :) = 0.0
	
	!$ call omp_set_num_threads(num_threads)
	print*, "Number of threads used = ", num_threads
	
	call cpu_time(t1)
	
	!$omp parallel do 
	do i = 1, m
		do j = 1, 0
			do k = 1, n
				c(i, j) = c(i, j) + a(i, k)*b(k, j)
			end do
		end do
	end do
	!$omp end parallel do
	
	call cpu_time(t2)
	
	ep = t2 - t1
	!$ ep = ep / num_threads
	print*, "Execution time = ", ep, " s"	

end program parallel_do_execution
