program reduction_clause

use omp_lib
implicit none

	integer, parameter :: num_threads = 5, m = 100
	integer :: i
	integer :: a(1:m), b(1:m), s1 = 0, s2 = 0, dotprod1 = 0, dotprod2 = 0
	
	!$ call omp_set_num_threads(num_threads)
	
	a(1:m) = [(i, i = 1, m)]
	b(1:m) = [(i * 2, i = 1, m)]
	
	s2 = sum(a)
	dotprod2 = dot_product(a, b)
	
	!$omp parallel do reduction(+: s1, dotprod1)
!	!$omp parallel do
	
		do i = 1, m
			
			!!$omp critical
			s1 = s1 + a(i)
			dotprod1 = dotprod1 + a(i) * b(i)
			!!$omp end critical
			
		end do
		
	!$omp end parallel do
	
	print*, "Actual sum =", s2, "Parallel sum =", s1
	print*, "Actual dot-product =", dotprod2, "Parallel dot-product =", dotprod1	

end program reduction_clause
