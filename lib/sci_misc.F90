module sci_misc

	use std_types
	use std_io

	interface gaussn
		module procedure gaussn_single, gaussn_multi
	end interface

contains

  !
  ! Calculates the expectation value of normally distributed x with mean m and standard deviation sd
  !
  subroutine gaussn_single(m,sd,x,val)            ! = Gaussian, normalised to unity
        real(dp), intent(in)  :: m
        real(dp), intent(in)  :: sd
        real(dp), intent(in)  :: x
        real(dp), intent(out) :: val

        val = (1.0_dp/(sd*sqrt(TWOPI_D)))*exp(-((x-m)**2.0_dp)/(2.0_dp*(sd**2.0_dp)))

  end subroutine gaussn_single

  !
  !
  !
  subroutine gaussn_multi(m,sd,x,val)            ! = Gaussian, normalised to unity
        real(dp), intent(in)  :: m
        real(dp), intent(in)  :: sd
        real(dp),dimension(:), intent(in)  :: x
        real(dp),dimension(size(x)), intent(out) :: val

        val = (1.0_dp/(sd*sqrt(TWOPI_D)))*exp(-((x-m)**2.0_dp)/(2.0_dp*(sd**2.0_dp)))

  end subroutine gaussn_multi

  !
  !  init random_real()
  !
  subroutine init_random_seed(s_e_e_d)
  		integer(i4b), intent(in), optional ::  s_e_e_d

		integer :: i, n, clock
		integer, dimension(:), allocatable :: seed

		call RANDOM_SEED(size = n)
		allocate(seed(n))

		call SYSTEM_CLOCK(COUNT=clock)

		seed = clock + 37 * (/ (i - 1, i = 1, n) /)
		if(present(s_e_e_d)) then
			seed = s_e_e_d
		end if
		call RANDOM_SEED(PUT = seed)

		deallocate(seed)
  end subroutine init_random_seed

  !
  !  return array of two times bigger size with the first array copied into the beginning
  !
  subroutine return_array_of_two_times_bigger_size(Ain, Aout)
 		real(dp), dimension(:), intent(in) :: Ain
 		real(dp), dimension(size(Ain)*2), intent(out) :: Aout

 		integer(i4b) :: i

 		Aout = 0.0_dp

 		do i=1,size(Ain)
 			Aout(i) = Ain(i)
 		end do
  end subroutine return_array_of_two_times_bigger_size

  !
  !   returns cross-product
  !
  function cross_product(a,b)
  		real(dp), dimension(:), intent(in) :: a,b
  		real(dp), dimension(size(a)) :: cross_product

  		if(size(a) /= size(b) .or. size(a) /= 3) then
  			call print_error_message(-1, "dimension error, arrays must have the same size of dimension 3")
  		end if


		cross_product(1) = -a(3)*b(2) + a(2)*b(3)
		cross_product(2) = a(3)*b(1) - a(1)*b(3)
		cross_product(3) = -a(2)*b(1) + a(1)*b(2)

  end function cross_product

  ! from http://orion.math.iastate.edu/burkardt/f_src/polpak/polpak.f90
  subroutine HermiteH ( n, x, cx )

    integer(i4b), intent(in)		:: n

    real(dp) cx(0:n)
    integer(i4b) i
    real(dp), intent(in)			::  x
!
    if ( n < 0 ) then
      return
    end if

    cx(0) = 1.0E+00

    if ( n == 0 ) then
      return
    end if

    cx(1) = 2.0E+00 * x

    do i = 2, n
      cx(i) = 2.0E+00 * x * cx(i-1) - 2.0E+00 * real ( i - 1 ) * cx(i-2)
    end do

    return
  end subroutine HermiteH

  function factorial (n) result (res)

    implicit none
    integer(i4b), intent (in) :: n
    real(dp) :: res
    integer(i4b) :: i

    res = product ((/(i, i = 1, n)/))

  end function factorial

  real(dp) recursive function franc_condon_factor(v,alpha,vv,aalpha,d) result(res)
  		integer(i4b), intent(in) 	:: v,vv
  		real(dp), intent(in)		:: alpha, aalpha, d

  		real(dp)					:: I, A, S, b, bb
  		integer(i4b)				:: k,kk,j

  		real(dp) herm1(0:v), herm2(0:vv)

		res = 0.0_dp

		S = alpha*aalpha*d*d/(alpha+aalpha)
		A = 2*sqrt(alpha*aalpha)/(alpha+aalpha)

		b = -aalpha*sqrt(alpha)*d/(alpha+aalpha)
		bb = alpha*sqrt(aalpha)*d/(alpha+aalpha)

		call HermiteH(v,b,herm1)
		call HermiteH(vv,bb,herm2)

  		do k=0,v
  		do kk=0,vv

  			if(mod(k+kk,2) == 1) then
  				cycle
  			end if

  			I = 1.0_dp
  			do j=1,k+kk-1,2
  				I = I*j
  			end do
  			I = I/((alpha+aalpha)**((k+kk)/2))

  			res = res + herm1(v-k)*herm2(vv-kk)*((2*sqrt(alpha))**k)*((2*sqrt(aalpha))**kk)*I*	&
  					factorial(v)/factorial(k)/factorial(v-k)*												&
  					factorial(vv)/factorial(kk)/factorial(vv-kk)

  		end do
  		end do

  		res = res * res * A * exp(-S)/(2**(v+vv))/factorial(v)/factorial(vv)
  		res = sqrt(res)

  end function franc_condon_factor


end module sci_misc
