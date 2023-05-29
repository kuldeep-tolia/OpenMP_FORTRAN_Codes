! OpenMP parallelized program for matrix multiplication of two square matrices using loop-unrolling method
program matrixMultiply

use omp_lib
implicit none

	integer :: num_thread = 1
	character(100) :: numchar
	character(100) :: nameq
	
	integer, parameter :: N = 3000
	integer :: i, j, k
	double precision, dimension(:, :), allocatable :: mat1, mat2, mat3
	double precision :: t1, t2
	
	if(COMMAND_ARGUMENT_COUNT() .ne. 1) then
		write(*, *) "Command line argument of number of threads is required."
		STOP
	end if
	
	call GET_COMMAND_ARGUMENT(0, nameq)	
	call GET_COMMAND_ARGUMENT(1, numchar)
	read(numchar, *) num_thread
	
	! Allocate matrices
	allocate(mat1(1:N, 1:N), mat2(1:N, 1:N), mat3(1:N, 1:N))
	
	! Initialize matrices
	mat1 = 1.0d0
	mat2 = 1.0d0
	mat3 = 0.0d0
	
	t1 = omp_get_wtime()
	!$omp parallel do num_threads(num_thread) default(shared) private(i, j, k)
	do j = 1, N, 2
		do k = 1, N
			do i = 1, N, 2		! Loop unrolling along with column access
				mat3(i, j) = mat3(i, j) + mat1(i, k) * mat2(k, j)
				mat3(i+1, j) = mat3(i+1, j) + mat1(i+1, k) * mat2(k, j)
				mat3(i, j+1) = mat3(i, j+1) + mat1(i, k) * mat2(k, j+1)
				mat3(i+1, j+1) = mat3(i+1, j+1) + mat1(i+1, k) * mat2(k, j+1)
			end do
		end do
	end do
	!$omp end parallel do
	t2 = omp_get_wtime()
	
	write(*, 10) "Number of threads =", num_thread
	10 format (A, I15)
	write(*, 10) "N =", N
	write(*, 11) "Execution time =", t2-t1, " seconds"
	11 format (A, F15.7, A)
	write(*, 12) "Sanity check ----> printing random index of the product matrix =", mat3(N/2, N/2)
	12 format (A, F15.3)
	
	! Deallocate matrices
	deallocate(mat1, mat2, mat3)
	
end program matrixMultiply
