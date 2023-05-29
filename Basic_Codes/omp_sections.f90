! OpenMP program to demonstrate the section construct
program omp_sections

use omp_lib
implicit none

	integer :: num_threads = 4, thread_id, i, j, k
	integer, parameter :: m = 10000, n = 10000
	real*8, dimension(1:m, 1:n) :: a, b, c, d, e
	
	b(:, :) = 0.0
	c(:, :) = 0.0
	d(:, :) = 0.0
	e(:, :) = 0.0
	do i = 1, m
		do j = 1, n
			a(i, j) = rand(1) * (3.0 * i + 0.1 * j)
		end do
	end do
	
	call omp_set_num_threads(num_threads)
	
	!$omp parallel
	!$omp sections private(i, j, thread_id)
	
		!$omp section
		thread_id = omp_get_thread_num()
		print*, "Section-1 started by thread id =", thread_id
		do j = 1, n
			do i = 1, m
				b(i, j) = a(i, j) * a(i, j)
			end do
		end do
		print*, "Section-1 ended by thread id =", thread_id
		
		!$omp section
		thread_id = omp_get_thread_num()
		print*, "Section-2 started by thread id =", thread_id
		do j = 1, n
			do i = 1, m
				c(i, j) = a(i, j) * 4.0 * (i + j)
			end do
		end do
		print*, "Section-2 ended by thread id =", thread_id
		
		!$omp section
		thread_id = omp_get_thread_num()
		print*, "Section-3 started by thread id =", thread_id
		do j = 1, n
			do i = 1, m
				d(i, j) = a(i, j)**2 + 1.0
			end do
		end do
		print*, "Section-3 ended by thread id =", thread_id
		
		!$omp section
		thread_id = omp_get_thread_num()
		print*, "Section-4 started by thread id =", thread_id
		do j = 1, n
			do i = 1, m
				e(i, j) = a(i, j) * 4.0
			end do
		end do
		print*, "Section-4 ended by thread id =", thread_id
		
	!$omp end sections
	!$omp end parallel

end program omp_sections
