!
! Standart types and constants of the NOSE package
!
!
module std_types

  integer, parameter :: i4b = selected_int_kind(9)
  integer, parameter :: i2b = selected_int_kind(4)
  integer, parameter :: i1b = selected_int_kind(2)
  integer, parameter :: sp  = kind(1.0)
  integer, parameter :: dp  = kind(1.0D0)
  integer, parameter :: spc = kind((1.0,1.0))
  integer, parameter :: dpc = kind((1.0D0,1.0D0))
  integer, parameter :: lgt = kind(.true.)
  real(sp), parameter :: PI=3.141592653589793238462643383279502884197_sp
  real(sp), parameter :: PIO2=1.57079632679489661923132169163975144209858_sp
  real(sp), parameter :: TWOPI=6.283185307179586476925286766559005768394_sp
  real(sp), parameter :: SQRT2=1.41421356237309504880168872420969807856967_sp
  real(sp), parameter :: EULER=0.5772156649015328606065120900824024310422_sp
  real(dp), parameter :: PI_D=3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: PIO2_D=1.57079632679489661923132169163975144209858_dp
  real(dp), parameter :: TWOPI_D=6.283185307179586476925286766559005768394_dp
  real(dp), parameter :: SQRT2_D=1.41421356237309504880168872420969807856967_dp
  real(dp), parameter :: EULER_D=0.5772156649015328606065120900824024310422_dp
  character(len=7), parameter  	:: FMT_D = 'E30.16'
  character(len=2), parameter	:: FMT_I = 'I8'

  !
  ! Here some important physical constants
  !

  real(dp), parameter :: kB_cmK = 0.6950387_dp

  ! energy conversion factor cm-1 -> internal
  real(dp), parameter :: Energy_cm_to_internal = 2.0_dp*PI_D*2.99792458d-5 ! 1.8836515673088532773...e-4
  real(dp), parameter :: Energy_internal_to_cm = 5308.8415534785_dp
  real(dp), parameter :: Energy_eV_to_interbal = 1.0_dp
  real(dp), parameter :: Energy_interval_to_eV = 1.0_dp !/Energy_eV_to_internal
  real(dp), parameter :: Energy_eV_to_cm = 8065.541_dp
  real(dp), parameter :: Energy_cm_to_eV = 1.0_dp/8065.541_dp
  real(dp), parameter :: kB_intK = kB_cmK * Energy_cm_to_internal ! 1.309210736595307880...e-4

  integer, dimension(:,:), allocatable, private :: std_type_kr


contains

	! initialization of the kroneker delta
	subroutine init_kroneker(N)
		integer :: N

		if (allocated(std_type_kr)) then
			deallocate(std_type_kr)
		end if

		allocate(std_type_kr(N,N))

		std_type_kr = 0.0d0
		do i = 1,N
			std_type_kr(i,i) = 1.0d0
		end do

	end subroutine init_kroneker

	! finalization of kroneker delta
	subroutine delete_kroneker

		if (allocated(std_type_kr)) then
			deallocate(std_type_kr)
		end if

	end subroutine delete_kroneker

	! kroneker delta
	function kroneker(i,j)
		integer :: i, j

		kroneker = std_type_kr(i,j)

	end function kroneker

end module std_types
