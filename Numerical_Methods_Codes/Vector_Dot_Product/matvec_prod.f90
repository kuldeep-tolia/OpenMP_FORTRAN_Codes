program mat_vec_prod

use omp_lib
implicit none

	integer, parameter :: n = 1000, m = 8
	double precision :: mat(n, n), b(n), c(n)
	double precision :: s, e, tmp
	integer :: i, j, k
	mat(:, :) = 1.0d0
	b(:) = 1.0d0
	call cpu_time(s)
	call omp_set_num_threads(m)
	!$omp parallel do default (shared) private (i, j, tmp)
	do i = 1, n
		tmp = 0.0d0
		do j = 1, n
			tmp = tmp + mat(j, i) * b(j)
		end do
		c(i) = tmp
	end do
	!$omp end parallel do
	
	call cpu_time(e)
	print*, "Number of threads used = ", m
	print*, "Time taken for operation = ", e-s, " s"
end program mat_vec_prod 
