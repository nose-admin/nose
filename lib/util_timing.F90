!
!
!
! TODO: Use std_io
!
module util_timing

	use std_types
    use std_io

	implicit none
	
	private 
	
	public :: init_timing, timing_clean
	public :: elapsed_time
	public :: next_code_position
	public :: print_characters
	public :: print_elapsed_time, print_total_elapsed_time
	
	public :: csec,cmin,chou,cday,cyea,proc_count,count_rate
	public :: count_max,code_position

	integer, dimension(:), allocatable  :: proc_count
	integer, dimension(:), allocatable  :: proc_count_store
	integer                    :: count_rate, count_max
	integer(i4b)               :: code_position
	integer(i4b)               :: position_max
	integer(i4b),parameter     :: position_max_incr = 100

	integer                  :: csec,cmin,chou,cday,cyea

contains

    
	subroutine init_timing()
		! initialize timing
		position_max = position_max_incr
		code_position = 1          ! start
		allocate(proc_count(position_max))
		call system_clock(proc_count(code_position),count_rate,count_max)
	end subroutine init_timing
	
	
	!
	! Advances the position in the code 
	!
	subroutine next_code_position()
		code_position = code_position + 1
		if (code_position > position_max) then
			! reallocate
			allocate(proc_count_store(position_max))
			proc_count_store = proc_count
			deallocate(proc_count)
			position_max = position_max + position_max_incr
			allocate(proc_count(position_max))
			proc_count(1:position_max-position_max_incr) = &
				proc_count_store
            deallocate(proc_count_store)
		end if
	end subroutine next_code_position
	
	!
	! Returns a string with human readable information about
	! the time elapsed between the last two code positions.
	!
	function elapsed_time(spec,cmax) result (str)
		character(len=*), intent(in)  :: spec
		logical, optional, intent(in) :: cmax
		character(len=256) :: str
		character(len=256) :: buffer
		integer(i4b) :: pos, ll
		real(dp)     :: rcsec
		
		!print *, " Code position: ",code_position		
		call system_clock(proc_count(code_position),count_rate,count_max)
		
		!print *, proc_count(code_position),proc_count(code_position-1)
		rcsec = real((proc_count(code_position)-proc_count(code_position-1))) &
			/count_rate
		csec = int(rcsec)
		if (present(cmax)) then
			if (cmax) then
				rcsec = real(count_max)/count_rate
				csec = int(rcsec)
			else
				rcsec = real((proc_count(code_position)- &
					proc_count(1)))/count_rate
				csec = int(rcsec)
			end if
		end if
	
		cmin = csec/60
		csec = csec - cmin*60
		rcsec = rcsec - cmin*60
		chou = cmin/60
		cmin = cmin - chou*60
		cday = chou/24
		chou = chou - cday*24
		cyea = cday/365
		cday = cday - cyea*365
	
		! here we have the number of years:  cyea
		!							  days:  cday
		!                            hours:  chou
		!                          minutes:  cmin
		!                      and seconds:  csec
		! that elapsed since last position in the code
		! we don't use years unless optional argument is specified and is true
				
		ll = len_trim(spec)
		str = ""

		do pos = 1, ll
		
			if (scan(spec,'s') == pos) then
				write(buffer,'(f8.4,a)') rcsec, " seconds "
			else if (scan(spec,'m') == pos) then
				write(buffer,'(i4,a)') cmin, " minutes "
			else if (scan(spec,'h') == pos) then
				write(buffer,'(i4,a)') chou, " hours "
			else if (scan(spec,'d') == pos) then
				write(buffer,'(i4,a)') cday, " days "
			else if (scan(spec,'y') == pos) then
				write(buffer,'(i4,a)') cyea, " years "
			else if (scan(spec,'a') == pos) then
				write(buffer,*) " and " 
			else
				buffer = ""
			end if
			str = trim(str)//trim(buffer)
			!print *, trim(str),trim(buffer)
		
		end do
		
	end function elapsed_time

	subroutine print_elapsed_time(chars)
		character(len=*), intent(in) :: chars
		
        call print_log_message("ELAPSED TIME - "//trim(chars),5) 

	end subroutine print_elapsed_time

	subroutine print_total_elapsed_time(chars)
		character(len=*), intent(in) :: chars
		
        call print_log_message("TOTAL ELAPSED TIME - "//trim(chars),5) 

	end subroutine print_total_elapsed_time

	subroutine print_characters(chars)
		character(len=*), intent(in) :: chars
		call print_log_message(trim(chars),5)
	end subroutine print_characters

	!
	!
	!
	subroutine timing_clean()
		deallocate(proc_count)
	end subroutine timing_clean



end module util_timing
