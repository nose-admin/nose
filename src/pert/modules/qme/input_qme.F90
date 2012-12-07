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
        real, dimension(:,:), allocatable :: rbuff
        integer, dimension(:,:), allocatable :: ibuff
        integer(i4b) :: rsize1, rsize2
        integer(i4b) :: isize1, isize2
        rsize1 = 100
        rsize2 = 100
        allocate(rbuff(rsize1,rsize2))
        isize1 = 100
        isize2 = 100
        allocate(ibuff(isize1,isize2))

        err = 0

        !
        ! excSource
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"input_qme")
            sz = nis_get_size()
            ibuff(1:4,1:1) = nis_get_integer(sz)
            Npos = ibuff(1:4,1)
        else
            call print_error_message(-1,"Missing record at NIS at position input_qme (excPosition)")
        end if
        wstr=adjustl("excSource        ")
        write(cbuff,'(a32,4i3)') wstr, Npos
        call print_log_message(trim(cbuff),5)
		
        !
        ! molecule
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"input_qme")
            sz = nis_get_size()
            write(cbuff,'(i3)') sz(1)
            write(molSystem,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
        else
            call print_error_message(-1,"Missing record at NIS at position input_qme (molecule)")
        end if
        wstr=adjustl("molecule         ")
        write(cbuff,'(a32,a)') wstr, trim(molSystem)
        call print_log_message(trim(cbuff),5)
		
        !
        ! pathway    
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"input_qme")
            sz = nis_get_size()
            write(cbuff,'(i3)') sz(1)
            write(pathway,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
        else
            call print_error_message(-1,"Missing record at NIS at position input_qme (pathway)")
        end if
        wstr=adjustl("pathway          ")
        write(cbuff,'(a32,a)') wstr, trim(pathway)
        call print_log_message(trim(cbuff),5)
		
        !
        ! number of runs (averaging over energy disorder)
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"input_qme")
            sz = nis_get_size()
            ibuff(1:2,1:1) = nis_get_integer(sz)
            Nruns   = ibuff(1,1)
            dgapini = ibuff(2,1)
        else
            call print_error_message(-1,"Missing record at NIS at position input_qme (runs)")
        end if
        wstr=adjustl("runs, dgapini    ")
        write(cbuff,'(a32,2i5)') wstr, Nruns, dgapini
        call print_log_message(trim(cbuff),5)
		
        !
        ! spec_fq (frequency boundary for the 2D spectrum)
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"input_qme")
            sz = nis_get_size()
            ibuff(1:3,1:1) = nis_get_integer(sz)
            spec_fq = ibuff(1:3,1)
        else
            call print_error_message(-1,"Missing record at NIS at position input_qme (spec_wini-dw-wst)")
        end if
        wstr=adjustl("spec_wini-dw-wst ")
        write(cbuff,'(a32,3i5)') wstr, spec_fq
        call print_log_message(trim(cbuff),5)
 		
        !
        ! localBasis
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"input_qme")
            sz = nis_get_size()
            write(cbuff,'(i3)') sz(1)
            write(locBasis,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
        else
            call print_error_message(-1,"Missing record at NIS at position input_qme (localBasis)")
        end if
        wstr=adjustl("localBasis       ")
        write(cbuff,'(a32,a)') wstr, trim(locBasis)
        call print_log_message(trim(cbuff),5)
		
        !
        ! relaxation
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"input_qme")
            sz = nis_get_size()
            write(cbuff,'(i3)') sz(1)
            write(doRelax,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
        else
            call print_error_message(-1,"Missing record at NIS at position input_qme (relaxation)")
        end if
        wstr=adjustl("relaxation       ")
        write(cbuff,'(a32,a)') wstr, trim(doRelax)
        call print_log_message(trim(cbuff),5)
		
        !
        ! secular approximation
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"input_qme")
            sz = nis_get_size()
            write(cbuff,'(i3)') sz(1)
            write(doSecAp,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
        else
            call print_error_message(-1,"Missing record at NIS at position input_qme (secular)")
        end if
        wstr=adjustl("secular ap.      ")
        write(cbuff,'(a32,a)') wstr, trim(doSecAp)
        call print_log_message(trim(cbuff),5)
		
        !
        ! pure dephasing calculation
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"input_qme")
            sz = nis_get_size()
            write(cbuff,'(i3)') sz(1)
            write(pureDephasing,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
        else
            call print_error_message(-1,"Missing record at NIS at position input_qme (dephasing)")
        end if
        wstr=adjustl("pure dephasing   ")
        write(cbuff,'(a32,a)') wstr, trim(pureDephasing)
        call print_log_message(trim(cbuff),5)
		
        !
        ! feeding rate
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"input_qme")
            sz = nis_get_size()
            rbuff(1:1,1:1) = nis_get_real(sz)
            Kfin = rbuff(1,1)
        else
            call print_error_message(-1,"Missing record at NIS at position input_qme (feeding)")
        end if
        wstr=adjustl("feeding          ")
        write(cbuff,'(a32,f9.5)') wstr, Kfin
        call print_log_message(trim(cbuff),5)
		
        !
        ! draining rate
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"input_qme")
            sz = nis_get_size()
            rbuff(1:1,1:1) = nis_get_real(sz)
            Kdin = rbuff(1,1)
        else
            call print_error_message(-1,"Missing record at NIS at position input_qme (draining)")
        end if
        wstr=adjustl("draining         ")
        write(cbuff,'(a32,f9.5)') wstr, Kdin
        call print_log_message(trim(cbuff),5)
		

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
