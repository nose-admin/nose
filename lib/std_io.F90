!
! Input/output routines
!
! Tomas Mancal  11/09/2005
!
!
module std_io

	use std_types

	implicit none

	!private

	public :: std_io_open
	public :: std_io_clean, init_std_io

    logical, public :: logging

	! interaces
	interface std_io_open
		module procedure open_file_1
	end interface

	! private variables
	logical, dimension(99) :: unit_pool, ext_open

	logical, private :: initiated

	integer(i4b), parameter :: reserved_units = 10
	integer(i4b), parameter :: max_units      = 99

    ! logging
    integer, private          :: loglevel
    character(len=9), private :: prompt
    character(len=9), private :: errprompt
    character(len=9), private :: warprompt

    ! output directory
    character(len=256)         :: out_dir

contains

	!
	! Initiates io and
    ! sets the log level
    !
    !
    ! Log levels: 0 .... no messages, not even errors
    !             5 .... standard level of NOSE
    !             10 ... all messages
    !
	subroutine init_std_io(p_id,log_level)
        integer, intent(in), optional :: p_id
        integer, intent(in), optional :: log_level
		unit_pool(1:reserved_units)  = .true.
		unit_pool(reserved_units+1:max_units) = .false.
		ext_open(1:reserved_units)   = .true.
		ext_open(reserved_units+1:max_units) = .false.
        if (present(log_level)) then
            loglevel  = log_level
        else
            loglevel  = 5
        end if
        prompt    = "Info   : "
        errprompt = "Error  : "
        warprompt = "Warning: "
        !if (parallel_id == 0) then
        if (present(p_id)) then
            if (p_id == 0) then
                logging = .true.
            else
                logging = .false.
            end if
        end if
  	end subroutine init_std_io

!
! FILE I/O
!

	!
	! Cleans before exiting (closes all open files)
	!
	subroutine std_io_clean()
		integer(i4b) :: i
		logical      :: op
		do i = reserved_units,max_units
			inquire(unit=i,opened=op)
			if (op.and.(.not.ext_open(i))) then
				 close(i)
				unit_pool(i) = .false.
			end if
		end do
	end subroutine

	!
	! opens file and returns unit number > 10
	!
	function open_file_1(filename) result (un)
		character(len=*), intent(in) :: filename
		integer                      :: un
		logical                      :: op

		integer(i4b) :: i

		un = 0
		if (.not.initiated) then
			call init_std_io()
		end if

		i = reserved_units
		do
			i = i + 1
			if (.not.unit_pool(i)) then
				inquire(unit=i,opened=op)
				if (op) then
					unit_pool(i) = .true.
					ext_open(i)  = .true.
				else
					un = i
					open(file=filename,unit=un)
					exit
				end if
			end if
		end do

		if (un == 0) stop 'No more free units'

	end function open_file_1

!
! LOGGING
!

    subroutine set_log_level(log_level)
        integer, intent(in) :: log_level
        loglevel = log_level
    end subroutine set_log_level

    function get_log_level() result (ll)
        integer :: ll
        ll = loglevel
    end function get_log_level

    !
    ! Priorities:
    !
    ! 1  ... highest
    !
    ! 10 ... lowest
    !

    !
    ! Prints message if its level is lower or equal to the log level
    !
    subroutine print_log_message(msg,ll)
        character(len=*), intent(in)  :: msg
        integer, intent(in)           :: ll
        if (ll <= loglevel) then
            if (.not.logging) return
            write(*,'(a,a)') prompt,msg
        end if

        call flush()
    end subroutine print_log_message

    !
    ! Prints message if its level is lower or equal to the log level
    !
    subroutine print_warning_message(msg,ll)
        character(len=*), intent(in)  :: msg
        integer, intent(in)           :: ll
        if (ll <= loglevel) then
            if (.not.logging) return
            write(*,'(a,a)') warprompt,msg
        end if

        call flush()
    end subroutine print_warning_message


    !
    ! Prints error
    !
    subroutine print_error_message(err,errmsg)
        integer, intent(in)           :: err
        character(len=*), intent(in)  :: errmsg
        if (err == 0) return
        if (loglevel > 0) then
            if (.not.logging) return
            write(*,'(a,a)') errprompt,errmsg
            stop
        end if

        call flush()
    end subroutine print_error_message

    !
    ! Returns true if this process can log with the given priority
    !
    function can_log_with(ll) result (r)
        logical :: r
        integer :: ll

        r = .false.
        if (ll <= loglevel) then
            if (.not.logging) return
            r = .true.
        end if

    end function can_log_with



    !
    ! Prints space
    !
    subroutine print_space(ll)
        integer, intent(in)           :: ll
        if (ll <= loglevel) then
            if (.not.logging) return
            write(*,'(a)') prompt
        end if
    end subroutine print_space


    !
    ! Dir name manipulation
    !
    function file_join(fdname1,fdname2) result (fdname3)
        character(len=*), intent(in) :: fdname1, fdname2
        character(len=256) :: fdname3

        if (len_trim(fdname1)>0) then
            fdname3 = trim(fdname1)//"/"//trim(fdname2)
        else
            fdname3 = trim(fdname2)
        end if

    end function file_join


    !
    ! Opens datafile and reads until a marker is reached
    !
    subroutine open_aux_file(unit,filename,marker,err)
    	integer :: unit
    	character(len=*) :: filename
    	character(len=*) :: marker

		character(len=128) :: line
		integer :: ios, err

    	open(unit=unit,file=filename,iostat=err)

    	if (err > 0) return

		err = 0

		do

			read(unit=unit,fmt=*,iostat=ios) line

			if (ios < 0) then
				err = 1   ! end of file reached
				print *, "End of file reached without finding the marker!"
				close(unit)
				exit
			end if

			if (trim(line).eq.trim(marker)) exit

		end do

    end subroutine open_aux_file

	subroutine close_aux_file(unit)
    	integer :: unit
		close(unit)

	end subroutine close_aux_file

end module std_io
