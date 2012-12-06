module input_montecarlo

    use resources_montecarlo

    implicit none

contains
    !
    ! Reads the input from the NOSE Input Stream (NIS)
    !
    subroutine read_input_montecarlo(err)
        integer, intent(out) :: err

        character(len=32) :: wstr
        character(len=256)::  cbuff
        real, dimension(:,:), allocatable :: rbuff
        integer, dimension(2) :: sz
        integer(i4b) :: i

        err = 0

        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"input_qme")
            sz = nis_get_size()
            write(cbuff,'(i3)') sz(1)
            write(methodMC,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
        else
            call print_error_message(-1,"Missing record at NIS at position input_qme")
        end if
        wstr=adjustl("moduleMethod     ")
        write(cbuff,'(a32,a)') wstr, trim(methodMC)
        call print_log_message(trim(cbuff),5)

        !
        ! Debug_gamma
        !
        allocate(rbuff(100,100))
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"reading NIS at position ? (debug gamma)")
            sz = nis_get_size()
            rbuff(1:1,1:1) = nis_get_real(sz)
            debug_gamma = rbuff(1,1)
        else
            call print_error_message(-1,"Missing record at NIS at position ? (debug gamma)")
        end if
        wstr=adjustl("debug gamma    ")
        write(cbuff,'(a32,f6.2)') wstr, debug_gamma
        call print_log_message(trim(cbuff),5)
        deallocate(rbuff)


		if(index(trim(methodMC), "_noG") == 0) then
			g_functions = .true.
		else
			g_functions = .false.
		end if

		if(index(trim(methodMC), "_fixed_seed") == 0) then
			fixed_seed = .false.
		else
			fixed_seed = .true.
		end if

		if(index(trim(methodMC), "_fastG") == 0) then
			vynechani_G_ifu = .false.
		else
			vynechani_G_ifu = .true.
		end if


        if (index(trim(methodMC), "normal") == 1) then

			depository = .false.
			variant360 = .false.
			modified_unraveling2 = .false.

        else if (index(trim(methodMC), "modified_unraveling") == 1) then

			depository = .false.
			variant360 = .false.
			modified_unraveling2 = .true.

        else if (index(trim(methodMC), "depository") == 1) then

			depository = .true.
			variant360 = .false.
			modified_unraveling2 = .false.

        else

			write(*,*) "Error: method not implemented, use normal[_noG], modified_unraveling[_noG], depository[_noG],"
			write(*,*) "       [_fixed_seed]"
            stop

        end if

    end subroutine read_input_montecarlo


end module input_montecarlo
