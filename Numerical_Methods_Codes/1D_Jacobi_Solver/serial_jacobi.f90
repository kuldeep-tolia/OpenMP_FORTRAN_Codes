program serial_jacobi

implicit none
	
	integer, parameter :: N = 10000		! Number of divisions
	double precision, parameter :: L = 1.0d0	! Length of domain 
	integer :: i, it_count
	double precision :: dx, total, comp_time, t1, t2
	double precision :: tolerance, global_error, residual 
	double precision, dimension(:), allocatable :: xgrid, phi, phi_new, phi_exact, b
	
	dx = L / N			! Grid numbering goes as 0:N
	tolerance = 1.0d-3		! Set the tolerance criteria
	it_count = 0
	residual = 0.0d0
	total = 0.0d0
	global_error = 1.0d0
	
	! Allocate vectors
	allocate(xgrid(0:N), phi(0:N), phi_new(0:N), phi_exact(0:N), b(0:N))
	
	! Initialize vectors
	do i = 0, N
		xgrid(i) = i * dx
		phi_exact(i) = ((-4.0d0 * xgrid(i)**3) + 1204.0d0 * xgrid(i) + 300.0d0) / 3.0d0
		b(i) = dx * dx * S(xgrid(i))
	end do
	
	phi = 100.0d0		! Initial guess value for scalar field
	phi_new = phi
	
	! Setting boundary conditions
	phi(0) = 100.0d0
	phi(N) = 500.0d0
	phi_new(0) = phi(0)
	phi_new(N) = phi(N)
	
	call cpu_time(t1)
	do while (global_error .ge. tolerance)
		
		total = 0.0d0
		
		! Jacobi iterative solver
		do i = 1, N-1
			phi_new(i) = (b(i) + phi(i-1) + phi(i+1)) / 2.0d0
		end do
		
		! Updating the vector
		do i = 1, N-1
			phi(i) = phi_new(i)
		end do
		
		! Check residual, r = b - Ax
		do i = 1, N-1
			residual = b(i) + phi(i-1) - 2.0d0 * phi(i) + phi(i+1)
			total = total + residual**2
		end do
		
		! Monitors
		global_error = sqrt(total)
		it_count = it_count + 1
		if (mod(it_count, 1000) .eq. 0) then
			write(*, 10) "Iteration =", it_count, " Global error =", global_error
			10 format (A, I15, A, F15.7)
		end if
	end do
	call cpu_time(t2)
	comp_time = t2 - t1
	
	! Post-processing
	open(11, file = "serial_solution_N_10000.txt")
	do i = 0, N
		write(11, *) xgrid(i), phi_exact(i), phi(i)
	end do
	
	write(*, *) "----------------------------"
	write(*, *) "N =", N
	write(*, 11) "Total time taken for convergence =", comp_time
	write(*, 11) "Final global error =", global_error
	11 format (A, F15.7)
	write(*, *) "Total iterations =", it_count

	deallocate(xgrid, phi, phi_new, phi_exact, b)
	
contains

	function S(x) result(y)
	implicit none
		
		double precision :: x, y
		y = 8.0d0 * x
		
	end function S

end program serial_jacobi
