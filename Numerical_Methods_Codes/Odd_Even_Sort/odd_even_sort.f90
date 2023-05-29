program tut3q2

use omp_lib
implicit none

	integer, parameter :: N1 = 10, N2 = 50
	integer :: i
	integer, dimension(:), allocatable :: vec1, vec2
	integer :: seed = 86234
	double precision :: r
	
	integer :: num_thread = 1
	character(100) :: numchar
	character(100) :: nameq
	
	if(COMMAND_ARGUMENT_COUNT() .ne. 1) then
		write(*, *) "Command line argument of number of threads is required."
		STOP
	end if
	
	call GET_COMMAND_ARGUMENT(0, nameq)	
	call GET_COMMAND_ARGUMENT(1, numchar)
	read(numchar, *) num_thread
	
	allocate(vec1(0:N1-1), vec2(0:N2-1))
	
	call random_seed(seed)
	do i = 0, N1-1
		call random_number(r)
		vec1(i) = floor(100 * r)
	end do
	
	do i = 0, N2-1
		call random_number(r)
		vec2(i) = floor(100 * r)
	end do
	
	write(*, *) "-----Unsorted vector-1-----"
	call displayVector(vec1, N1)
	call oddEvenSort(vec1, N1)
	write(*, *) "-----Sorted vector-1-----"
	call displayVector(vec1, N1)
	
	write(*, *) "-----Unsorted vector-2-----"
	call displayVector(vec2, N2)
	call oddEvenSort(vec2, N2)
	write(*, *) "-----Sorted vector-2-----"
	call displayVector(vec2, N2)	
		
	deallocate(vec1, vec2)
	
contains

	subroutine displayVector(vec, n)
	
	implicit none
	
		integer :: n, i
		integer, dimension(0:n-1) :: vec
		
		do i = 0, n-1
			write(*, *) vec(i)
		end do
		
		write(*, *) " "
	
	end subroutine displayVector
	
	subroutine oddEvenSort(vec, n)
	
	implicit none
	
		integer :: i, n, phase
		integer :: thread_count
		integer, dimension(0:n-1) :: vec
		
		thread_count = omp_get_num_threads()
		
		!$omp parallel num_threads(thread_count) default(none) shared(vec, n) private(i, phase)
		do phase = 0, n-1
			if (mod(phase, 2) .eq. 0) then
				!$omp do
				do i = 1, n-1, 2
					if (vec(i-1) .gt. vec(i)) then
						vec(i) = vec(i) + vec(i-1)
						vec(i-1) = vec(i) - vec(i-1)
						vec(i) = vec(i) - vec(i-1)
					end if
				end do
				!$omp end do
			else
				!$omp do
				do i = 1, n-2, 2
					if (vec(i) .gt. vec(i+1)) then
						vec(i) = vec(i) + vec(i+1)
						vec(i+1) = vec(i) - vec(i+1)
						vec(i) = vec(i) - vec(i+1)
					end if
				end do
			end if
		end do
		!$omp end parallel
	
	end subroutine oddEvenSort	
	
end program tut3q2
