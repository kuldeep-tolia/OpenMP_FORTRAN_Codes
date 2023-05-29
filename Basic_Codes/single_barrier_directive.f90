program single_barrier_directive

use omp_lib
implicit none

	integer :: num_threads = 16, istart, iend, thread_id, i, points_per_thread
	integer, parameter :: m = 1000000000
	double precision, parameter :: pi = 4.0d0 * atan(1.0d0)
	double precision :: t1, t2, ep, sum1, sum2, x, dx, l = 90.0d0
	
	dx = l * pi / (180.0d0 * (m - 1))
	points_per_thread = (m + num_threads - 1.0d0) / float(num_threads)
	istart = 1
	iend = m
	sum1 = 0.0d0
	sum2 = 0.0d0
	
	!$ call omp_set_num_threads(num_threads)
	
	call cpu_time(t1)
	
	!$omp parallel default(shared) private(i, thread_id, x, istart, iend, sum1)
	
		!$ thread_id = omp_get_thread_num()
		!$ istart = thread_id * points_per_thread + 1
		!$ iend = min(m, thread_id * points_per_thread + points_per_thread)
		
		do i = istart, iend
			x = (i - 1) * dx
			sum1 = sum1 + sin(x) * dx
		end do
		
		!$omp barrier
		
		!$omp critical
			sum2 = sum2 + sum1
		!$omp end critical
		
		!$omp single
			!$ print*, "Single construct performed by thread id = ", thread_id
			print*, "Area of integral = ", sum2
		!$omp end single
		
	!$omp end parallel
	
	call cpu_time(t2)
	ep = t2 - t1
	!$ ep = ep / num_threads
	
	!$ print*, "OpenMP mode ON!!"
	print*, "Time taken for computation = ", ep, " s"	

end program single_barrier_directive
