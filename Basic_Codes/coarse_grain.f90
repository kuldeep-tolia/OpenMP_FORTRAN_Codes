program coarse_grain_parallelism

use omp_lib
implicit none

	integer :: num_threads = 16, istart, iend, thread_id, i, points_per_thread
	integer, parameter :: m = 1000000000
	double precision, dimension(1:m) :: a
	double precision, parameter :: pi = 4.0d0 * atan(1.0d0)
	double precision :: t1, t2, elapsed_time, x, dx, l = 10.0d0
	
	a(:) = 0.0d0		!Initialisation
	dx = 1.0d0 / float(m - 1)
	points_per_thread = (m + num_threads - 1.0d0) / float(num_threads)
	istart = 1
	iend = m
	
	!$ call omp_set_num_threads(num_threads)
	t1 = omp_get_wtime()
	
	!$omp parallel default(shared) private(istart, iend, i, x, thread_id)
	
		thread_id = omp_get_thread_num()
		istart = thread_id * points_per_thread + 1
		iend = min(m, thread_id * points_per_thread + points_per_thread)
		
		do i = istart, iend
			x = (i - 1) * dx
			a(i) = (sin(2 * pi * x / l))**2 + (cos(2 * pi * x / l))**2
		end do
		
		!$omp critical
			print*, "thread id = ", thread_id, " istart = ", istart, " iend = ", iend
		!$omp end critical
	
	!$omp end parallel
	
	t2 = omp_get_wtime()
	elapsed_time = t2 - t1
	
	print*, "OpenMP Mode ON"
	print*, "Values of 'a' are = ", a
	print*, "Time taken for computation = ", elapsed_time, " s"

end program coarse_grain_parallelism
