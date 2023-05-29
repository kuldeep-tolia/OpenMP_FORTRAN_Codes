! OpenMP program to solve given Poisson equation using Gausss-Seidel method with diagonal way of traversing grid points
program openmp_diagonal_gauss_seidel

use omp_lib
implicit none
	
	integer :: num_thread = 1
	character(100) :: numchar
	character(100) :: nameq
	
	double precision, parameter :: delta = 1.0d-2	! grid spacing, assumed delta_x = delta_y = delta
	integer :: i, j, l, N, istart, iend, it_count
	double precision :: xmin, xmax, ymin, ymax
	double precision :: tolerance, global_error, t1, t2
	double precision, dimension(:), allocatable :: xgrid, ygrid
	double precision, dimension(:, :), allocatable :: phi, q, phi_exact, domain_error
	
	! Set domain extent
	xmin = -1.0d0
	xmax = 1.0d0
	ymin = -1.0d0
	ymax = 1.0d0
	
	N = (xmax - xmin) / delta		! Number of grid points required in x, y directions	NOTE: grid numbering goes as 0:N
	tolerance = 1.0d-6		! Set tolerance value
	global_error = 1.0d0
	it_count = 0
	
	if(COMMAND_ARGUMENT_COUNT() .ne. 1) then
		write(*, *) "Command line argument of number of threads is required."
		STOP
	end if
	
	call GET_COMMAND_ARGUMENT(0, nameq)	
	call GET_COMMAND_ARGUMENT(1, numchar)
	read(numchar, *) num_thread
	
	write(*, 11) "N =", N, " delta =", delta
	11 format (A, I7, A, F15.5)
	write(*, *) "Number of threads =", num_thread
	
	! Allocate vectors and matrices
	call allocate_memory()
	
	t1 = omp_get_wtime()
	! Fill in grid locations for x, y discrete points
	call initialize_grid()		
	
	! Calculate phi_exact with the provided solution
	call calc_phi_exact()
	
	! Fill in source term (q) 
	call calc_source_term()
	
	! Initialize the domain with guess value
	phi(:, :) = 0.0d0
	
	! Set boundary conditions
	do i = 0, N
		phi(i, 0) = 0.0d0		! Bottom boundary
		phi(i, N) = 0.0d0		! Top boundary
	end do
	
	do j = 0, N
		phi(0, j) = 0.0d0		! Left boundary
		phi(N, j) = 0.0d0		! Right boundary
	end do
	
	!$omp parallel num_threads(num_thread) default(shared) private(i, j)
	do while (global_error .ge. tolerance)
		
		! Gauss-Seidel iterative solver with diagonal way of traversing grid points
		do l = 0, 2*N
			if (l .le. N) then
				istart = 0
				iend = l
			else
				istart = l - N
				iend = N
			end if
			
			!$omp do
			do i = istart, iend
				j = l - i
				if ((i .gt. 0) .and. (i .lt. N) .and. (j .gt. 0) .and. (j .lt. N)) then
					phi(i, j) = 0.250d0 * (phi(i+1, j) + phi(i-1, j) + phi(i, j+1) + &
					&	   phi(i, j-1)) + (q(i, j) * 0.250d0 * delta**2)
					
					domain_error(i, j) = abs((phi_exact(i, j) - phi(i, j)) / &
					&		   phi_exact(i, j))
				end if
			end do
			!$omp end do
		end do
		
		! Calculate global error by taking max-norm of error-matrix
		!$omp single
		global_error = maxval(domain_error)
		
		! Monitors
		it_count = it_count + 1
		if (mod(it_count, 10) .eq. 0) then
			write(*, 10) "#Iteration =", it_count, " Global error =", global_error
			10 format (A, I15, A, F15.7)
		end if
		!$omp end single 
	end do
	!$omp end parallel
	t2 = omp_get_wtime()
	
	! Post-processing
	open(11, file = "phi_profile_y_05.txt")
	do i = 0, N
		write(11, *) xgrid(i), phi_exact(i, 3*N/4), phi(i, 3*N/4)
	end do
	
	write(*, *) "----------------------------"
	write(*, 12) "Total time taken for convergence =", t2-t1
	write(*, 12) "Final global error =", global_error
	12 format (A, F15.7)
	write(*, *) "Total iterations taken for convergence =", it_count
	
	call deallocate_memory()
	
contains
	subroutine allocate_memory()
	implicit none
		
		allocate (xgrid(0:N), ygrid(0:N))
		allocate (phi(0:N, 0:N), q(0:N, 0:N), phi_exact(0:N, 0:N), domain_error(0:N, 0:N))
		
	end subroutine allocate_memory
	
	subroutine initialize_grid()
	implicit none
		
		xgrid(0) = xmin
		ygrid(0) = ymin
		do i = 1, N
			xgrid(i) = xgrid(i-1) + delta
			ygrid(i) = ygrid(i-1) + delta
		end do
		
	end subroutine initialize_grid
	
	subroutine calc_phi_exact()
	implicit none
		
		!$omp parallel do if(N .ge. 100) num_threads(num_thread) default(shared) private(i, j) collapse(2)
		do j = 0, N
			do i = 0, N
				phi_exact(i, j) = (xgrid(i)**2 - 1.0d0) * (ygrid(j)**2 - 1.0d0)
			end do
		end do
		!$omp end parallel do
	
	end subroutine calc_phi_exact
	
	subroutine calc_source_term()
	implicit none
		
		!$omp parallel do if(N .ge. 100) num_threads(num_thread) default(shared) private(i, j) collapse(2)
		do j = 0, N
			do i = 0, N
				q(i, j) = 2.0d0 * (2.0d0 - xgrid(i)**2 - ygrid(j)**2)
			end do
		end do
		!$omp end parallel do
		
	end subroutine calc_source_term
	
	subroutine deallocate_memory()
	implicit none
		
		deallocate (xgrid, ygrid)
		deallocate (phi, q, phi_exact, domain_error)
		
	end subroutine deallocate_memory

end program openmp_diagonal_gauss_seidel
