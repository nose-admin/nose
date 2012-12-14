!
!
!
module util_allocation

  use std_types
  use iso_c_binding

  !integer, private, parameter :: c_size_t = i4b
  integer(i8b), private :: current_byte_size
  integer(i8b), private :: total_byte_size
  integer(i8b), private :: max_byte_size

  interface util_sizeof
  	module procedure util_sizeof_i,util_sizeof_ii
  end interface

contains

  subroutine init_allocation_counting()
    current_byte_size = 0
    total_byte_size   = 0.0_dp
    max_byte_size     = 0.0_dp
  end subroutine init_allocation_counting

  subroutine clean_allocation_counting()
    call init_allocation_counting()
  end subroutine clean_allocation_counting

  function get_byte_count() result(r)
    integer(i8b) :: r
    r = total_byte_size
  end function get_byte_count

  function get_mb_count() result(r)
    real(dp) :: r
    r = real(total_byte_size,dp)/(1024.0_dp**2)
  end function get_mb_count

  function get_max_byte_count() result(r)
    integer(i8b) :: r
    r = max_byte_size
  end function get_max_byte_count

  function get_max_mb_count() result(r)
    real(dp) :: r
    r = real(max_byte_size,dp)/(1024.0_dp**2)
  end function get_max_mb_count

  subroutine count_bytes(arg)
    integer(c_size_t) :: arg
    total_byte_size = total_byte_size + arg
    if (total_byte_size > max_byte_size) then
      max_byte_size = total_byte_size
    end if
    !print *, real(arg)/real(1024*1024)
    print *, int(real(get_byte_count())/real(1024**2)),'MB'
    call flush()
  end subroutine count_bytes

  subroutine uncount_bytes(arg)
    integer(c_size_t) :: arg
    total_byte_size = total_byte_size - arg
    print *, int(real(get_byte_count())/real(1024**2)),'MB'
    call flush()
  end subroutine uncount_bytes

  subroutine count_name(name)
    character(len=*) :: name
    write(*,'(A,A,I6,A)', advance='no') "Allocating ", name, int(real(get_byte_count())/real(1024**2)),' MB  ->'
    call flush()
  end subroutine count_name

  subroutine uncount_name(name)
    character(len=*) :: name
	write(*,'(A,A,I6,A)', advance='no') "Deallocating ", name, int(real(get_byte_count())/real(1024**2)),'MB'
    call flush()
  end subroutine uncount_name

  function util_sizeof_i(i) result (res)
  	integer :: i
  	integer(i8b) :: res
  	res = 2*i
  end function util_sizeof_i

  function util_sizeof_ii(i,j) result (res)
  	integer :: i,j
  	integer(i8b) :: res
  	res = util_sizeof_i(i)*j
  end function util_sizeof_ii

end module util_allocation
