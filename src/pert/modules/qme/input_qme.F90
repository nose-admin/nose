module input_qme

    use resources_qme

    implicit none

contains
    !
    ! Reads the input from the NOSE Input Stream (NIS)
    !
    subroutine read_input_qme(err)
        integer, intent(out) :: err

        character(len=32) :: wstr
        character(len=256)::  cbuff
        integer, dimension(2) :: sz
        integer(i4b) :: i

        err = 0

        !
        ! use twoexc.
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"input_qme")
            sz = nis_get_size()
            write(cbuff,'(i3)') sz(1)
            write(method,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
        else
            call print_error_message(-1,"Missing record at NIS at position input_qme")
        end if
        wstr=adjustl("moduleMethod     ")
        write(cbuff,'(a32,a)') wstr, trim(method)
        call print_log_message(trim(cbuff),5)



		if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"input_qme")
            Nte = nis_get_integer()

            ! resolve consequences
            grid_Nt = max(grid_Nt,Nte*gt(1))
			deallocate(grid_t)
			allocate(grid_t(grid_Nt))                                 ! MARK for deallocation
			do i = 1, grid_Nt
            	grid_t(i) = real((i-1),dp)*dt
			end do

        else
            call print_error_message(-1,"Missing record at NIS at position input_qme")
        end if











		if(index(trim(method), "_normU") == 0) then
			tau_projector_normalization_for_others = .false.
		else
			tau_projector_normalization_for_others = .true.
		end if

        if (index(trim(method), "PT2-SBsm") == 1 .or. index(trim(method), "PT2-SBssm") == 1) then

            use_twoexcitons = .true.
            read_external_evops = .false.
            use_module_nakajima_zwanzig = .true.
            submethod1 = 's'
            submethod2 = 's'

        else if (index(trim(method), "PT2-SBssM") == 1) then

            use_twoexcitons = .true.
            read_external_evops = .false.
            use_module_nakajima_zwanzig = .true.
            submethod1 = 's'
            submethod2 = 'S'

        else if (index(trim(method), "PT2-SBsM") == 1) then

            use_twoexcitons = .true.
            read_external_evops = .false.
            use_module_nakajima_zwanzig = .true.
            submethod1 = 'S'
            submethod2 = 'S'

        else if (index(trim(method), "PT2-SBsq") == 1) then

            use_twoexcitons = .true.
            read_external_evops = .false.
            use_module_nakajima_zwanzig = .true.
            submethod1 = 'Q'
            submethod2 = 'Q'

        else if (index(trim(method), "PT2-SBfm") == 1) then

            use_twoexcitons = .true.
            read_external_evops = .false.
            use_module_nakajima_zwanzig = .true.
            submethod1 = 'r'
            submethod2 = 'r'

        else if (index(trim(method), "PT2-SBsfm") == 1) then

            use_twoexcitons = .true.
            read_external_evops = .false.
            use_module_nakajima_zwanzig = .true.
            submethod1 = 's'
            submethod2 = 'r'

        else if (index(trim(method), "PT2-SBsmp") == 1) then ! 's' + phenomenological relaxation in E block

            use_twoexcitons = .true.
            read_external_evops = .false.
            use_module_nakajima_zwanzig = .true.
            submethod1 = 'r'
            submethod2 = 'p'

        else if (index(trim(method), "PT2-SBsMp") == 1) then ! 'S' + phenomenological relaxation in E block

            use_twoexcitons = .true.
            read_external_evops = .false.
            use_module_nakajima_zwanzig = .true.
            submethod1 = 'r'
            submethod2 = 'P'

        else if (index(trim(method), "PT2-SBfM") == 1) then

            use_twoexcitons = .true.
            read_external_evops = .false.
            use_module_nakajima_zwanzig = .true.
            submethod1 = 'R'
            submethod2 = 'R'

        else if (index(trim(method), "PT2-SBsfM") == 1) then

            use_twoexcitons = .true.
            read_external_evops = .false.
            use_module_nakajima_zwanzig = .true.
            submethod1 = 's'
            submethod2 = 'R'

        else if (index(trim(method), "PT2-SBfq") == 1) then

            use_twoexcitons = .true.
            read_external_evops = .false.
            use_module_nakajima_zwanzig = .true.
            submethod1 = 'q'
            submethod2 = 'q'

        else if (index(trim(method), "PT2-SBsfq") == 1) then

            use_twoexcitons = .true.
            read_external_evops = .false.
            use_module_nakajima_zwanzig = .true.
            submethod1 = 's'
            submethod2 = 'q'

        else if (index(trim(method), "PT2-SBssq") == 1) then

            use_twoexcitons = .true.
            read_external_evops = .false.
            use_module_nakajima_zwanzig = .true.
            submethod1 = 's'
            submethod2 = 'Q'

        else if (index(trim(method), "PT2-SBsfu") == 1) then

            use_twoexcitons = .true.
            read_external_evops = .false.
            use_module_nakajima_zwanzig = .true.
            submethod1 = 's'
            submethod2 = 'u'

        else if (index(trim(method), "PT2-SBssu") == 1) then

            use_twoexcitons = .true.
            read_external_evops = .false.
            use_module_nakajima_zwanzig = .true.
            submethod1 = 's'
            submethod2 = 'U'

        else if (index(trim(method), "PT2-SBsSPEC") == 1) then

            use_twoexcitons = .true.
            read_external_evops = .false.
            use_module_nakajima_zwanzig = .true.
            submethod1 = 's'
            submethod2 = '#'

        else if (index(trim(method), "carotenoid-special") == 1) then

            use_twoexcitons = .true.
            read_external_evops = .false.
            use_module_nakajima_zwanzig = .true.
            submethod1 = '-'
            submethod2 = '-'

        else if (trim(method) == "PT2-RC") then

            use_twoexcitons = .false.
            read_external_evops = .false.
            use_module_nakajima_zwanzig = .false.

        else

			write(*,*) "Error: method not implemented, use PT2-RC, PT2-SBs[sf][mMqu] or PT2-SB[sf][mMqu]."
            stop

        end if


    end subroutine read_input_qme


end module input_qme
