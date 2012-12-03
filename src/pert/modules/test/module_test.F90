module module_test
    
	use prepare
	use twod

	use std_types

	use qch_lib

	use numer_ode
	use numer_fft
	use numer_matrix
	use sci_misc

	implicit none

	contains

	subroutine do_test_work(err)
		integer, intent(out) :: err

		real(dp), dimension(2,2) 		:: S, H, V, T
		real(dp), dimension(2,2,2,2) 	:: four_center
		real(dp)						:: E
		integer(i4b)					:: i, j, k, l

		call test_hydrogen_molecule(four_center, S, H, T, V, E)

		err = 0

		open(unit=11,file=trim(file_join(out_dir,"hydrogen_molecule.dat")))

		do i=1,2
		do j=1,2
			write(11,*) S(i,j)
		end do
		end do
		do i=1,2
		do j=1,2
			write(11,*) H(i,j)
		end do
		end do
		do i=1,2
		do j=1,2
			write(11,*) T(i,j)
		end do
		end do
		do i=1,2
		do j=1,2
			write(11,*) V(i,j)
		end do
		end do
		do i=1,2
		do j=1,2
		do k=1,2
		do l=1,2
			write(11,*) four_center(i,j,k,l)
		end do
		end do
		end do
		end do

		write(11,*) E

		close(11)


	end subroutine do_test_work

	subroutine collect_test_data(err)
		integer, intent(out) :: err

		err = 0


	end subroutine collect_test_data

end module module_test
