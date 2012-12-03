!
! Module for printing nice formated tables of data into standard output
!
! Tomas Mancal
!
module util_char_graphics

	use std_types
    use std_io

	!private

	! types
	public :: table

	! constructor
	public :: new_table

	! destructor
	public :: delete_table

	! methods
	public :: next_tableLine


	!*******************************

	!
	! Public types
	!
	type table
		integer :: tab_nr        ! index of the table in the internal array
	end type table

	!
	! Subroutine interfaces
	!
	interface next_tableLine
		module procedure next_tableLine_dp
	end interface

	!****************
	! Internals
	!****************

	! basic table
	type i_table
		integer(i4b) :: N_columns
		integer(i4b) :: N_lines

	end type i_table

	! holds the current table
	type(i_table), pointer :: current_table

	! holds all tables
	type(i_table), dimension(:), allocatable, target :: tables

	! current number of tables
	integer(i4b) :: N_tables

	! actual number of tables hold
	integer(i4b) :: N_tab_actual

contains


    subroutine print_devider(i)
        integer, intent(in) :: i

        if (i == 1) then
            call print_log_message(&
"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++",5)
        else if (i == 2) then
            call print_log_message(&
"======================================================================",5)
        else if (i == 3) then
            call print_log_message(&
"----------------------------------------------------------------------",5)
        end if

    end subroutine print_devider

	!
	! creates new table object
	!
	function new_table() result (this)
		type(table) :: this
		this%tab_nr = 1
		if (N_tables == 0) then
			N_tables = 1
			N_tab_actual = 1
			allocate(tables(N_tab_actual))
			tables(1)%N_columns = 10
			tables(2)%N_lines   = 1
			current_table => tables(this%tab_nr)
		end if

		print *, current_table%N_columns

	end function

	!
	! destroy the table object
	!
	subroutine delete_table(this)
		type(table), intent(inout) :: this
		this%tab_nr = -1
	end subroutine


	!**********
	! Methods
	!**********


	subroutine next_tableLine_dp(vals)
		real(dp), dimension(:), intent(in) :: vals

	end subroutine next_tableLine_dp

end module util_char_graphics

