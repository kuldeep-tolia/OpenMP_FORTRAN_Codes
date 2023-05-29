program tut4q3

use omp_lib
implicit none

	integer :: num_thread = 1
	character(100) :: numchar
	character(100) :: nameq
	
	integer :: ndiv, n, j
	double precision :: xstart, xend, dx
	double precision, dimension(:), allocatable :: u, xgrid, dudx_exact
	double precision, dimension(:, :), allocatable :: dudx
	
	dx = 0.0010d0		! Set delta-x parameter
	xstart = -1.0d0		! Lower limit of domain
	xend = 1.0d0		! Upper limit of domain
	
	ndiv = (xend - xstart) / dx		! Number of divisions
	n = ndiv+1				! Number of grid points (1:n)
	
	if(COMMAND_ARGUMENT_COUNT() .ne. 1) then
		write(*, *) "Command line argument of number of threads is required."
		STOP
	end if
	
	call GET_COMMAND_ARGUMENT(0, nameq)	
	call GET_COMMAND_ARGUMENT(1, numchar)
	read(numchar, *) num_thread
	
	! Allocate vectors
	allocate(u(1:n), xgrid(1:n), dudx_exact(1:n))
	allocate(dudx(1:n, 1:2))
	
	call initialize_vectors(n, xgrid, u, dudx_exact, dudx)	! Initialize vectors
	call cds_2order(n, dudx(:, 1), u, num_thread)
	call cds_4order(n, dudx(:, 2), u, num_thread)
	
	! Post-processing
	open(11, file = "cds_dx2.txt")
	do j = 1, n
		write(11, *) xgrid(j), dudx_exact(j), dudx(j, 1), dudx(j, 2)
	end do
	
	deallocate(dudx)
	deallocate(u, xgrid, dudx_exact)
	
contains

	function f(x) result(y)
	
	implicit none
		
		double precision :: x, y
		
		y = 7.0d0 - x * tan(x)
		
	end function f
	
	function df(x) result(y)
	
	implicit none
		
		double precision :: x, y
		
		y = -tan(x) - x / (cos(x)**2)
		
	end function df
	
	subroutine initialize_vectors(m, v1, v2, v3, v4)
	
	implicit none
	
		integer :: m, j
		double precision, dimension(1:m) :: v1, v2, v3
		double precision, dimension(1:m, 1:2) :: v4
		
		do j = 1, m
			v1(j) = xstart + (j-1) * dx
			v2(j) = f(v1(j))
			v3(j) = df(v1(j))
			v4(j, 1) = 0.0d0
			v4(j, 2) = 0.0d0
		end do
	
	end subroutine initialize_vectors
	
	subroutine cds_2order(m, dzdx, u, nthread)
	
	implicit none
	
		integer :: m, nthread, j
		double precision, dimension(1:m) :: dzdx, u
		
		! Calculate first derivative using forward/backward difference formulae near boundary points
		dzdx(1) = (u(2) - u(1)) / dx	
		dzdx(m) = (u(m) - u(m-1)) / dx
		
		! Calculate first derivative using central difference formula 2nd order for interior points
		!$omp parallel do num_threads(nthread) default(shared) private(j)
		do j = 2, m-1
			dzdx(j) = (u(j+1) - u(j-1)) / (2.0d0 * dx)
		end do
		!$omp end parallel do
	
	end subroutine cds_2order
	
	subroutine cds_4order(m, dzdx, u, nthread)
	
	implicit none
	
		integer :: m, nthread, j
		double precision, dimension(1:m) :: dzdx, u
		
		! Calculate first derivative using forward/backward/central difference formulae near boundary points
		dzdx(1) = (u(2) - u(1)) / dx;
		dzdx(2) = (u(3) - u(1)) / (2.0d0 * dx)	
		dzdx(m) = (u(m) - u(m-1)) / dx
		dzdx(m-1) = (u(m) - u(m-2)) / (2.0d0 * dx)
		
		! Calculate first derivative using central difference formula 4th order for interior points
		!$omp parallel do num_threads(nthread) default(shared) private(j)
		do j = 3, m-2
			dzdx(j) = (8.0d0 * u(j+1) - 8.0d0 * u(j-1) + u(j-2) - u(j+2)) / (12.0d0 * dx)
		end do
		!$omp end parallel do
	
	end subroutine cds_4order

end program tut4q3
