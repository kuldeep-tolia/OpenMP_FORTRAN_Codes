! OpenMP program to demonstrate the master construct
program master_construct

use omp_lib
implicit none

	integer, parameter :: num_threads = 4, n = 12
	integer :: i, istart, iend, thread_id, points_per_thread
	real*8 :: a = 50.0, b(1:n) = 0.0
	
	points_per_thread = (n + num_threads - 1) / num_threads
	istart = 1
	iend = n
	
	call omp_set_num_threads(num_threads)
	
	!$omp parallel default(shared) private(i, istart, iend, thread_id)
		
		!$omp master
			!$ thread_id = omp_get_thread_num()
			!$ a = real(thread_id) + 100.0
			!$ print*, "Master block working here by thread number =", thread_id
		!$omp end master
		!$omp barrier
		thread_id = omp_get_thread_num()
		istart = thread_id * points_per_thread + 1
		iend = min(n, thread_id * points_per_thread + points_per_thread)
		print*, "thread number =", thread_id, "istart =", istart, "iend =", iend
		
		do i = istart, iend
			b(i) = a
		end do
		
	!$omp end parallel
	
	do i = 1, n
		print '("Array-element b(",i2,") = ",f15.6)', i, b(i)
	end do

end program master_construct
