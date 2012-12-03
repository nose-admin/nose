!
! Reads and contains all information from simulate.in file
!
module simulate_wizard

    use resources

    ! module in use
    use resources_tdpt3
    use resources_qme

    implicit none

!    private :: module_name

!    character(len=64) :: module_name

    character(len=256), private :: outp
    character(len=32), private  :: wstr

contains

    !
    ! Initialization and reading of the simulation parameters
    !
    subroutine init_simulation(err)
        integer, intent(out) :: err
        real, dimension(:,:), allocatable :: rbuff
        integer, dimension(:,:), allocatable :: ibuff
        character(len=256) :: cbuff
        integer(i4b) :: rsize1, rsize2
        integer(i4b) :: isize1, isize2

        integer, dimension(2) :: sz

        integer :: i

        err = 0

        rsize1 = 100
        rsize2 = 100
        allocate(rbuff(rsize1,rsize2))
        isize1 = 100
        isize2 = 100
        allocate(ibuff(isize1,isize2))

        ! check if NIS is ready and start reading from it
        call nis_rewind()

        call print_log_message(" ",5)
        call print_log_message("--- Configuration ---",5)
        !
        ! Module name
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"reading NIS at position 1 (module name)")
            sz = nis_get_size()
            write(cbuff,'(i3)') sz(1)
            write(modname,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
        else
            call print_error_message(-1,"Missing record at NIS at position 1 (module name)")
        end if
        wstr=adjustl("moduleName     ")
        write(cbuff,'(a32,a)') wstr, trim(modname)
        call print_log_message(trim(cbuff),5)

        !
        ! Temperature
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"reading NIS at position 2 (temperature)")
            sz = nis_get_size()
            rbuff(1:1,1:1) = nis_get_real(sz)
            temp = rbuff(1,1)
        else
            call print_error_message(-1,"Missing record at NIS at position 2 (temperature)")
        end if
        wstr=adjustl("temperature    ")
        write(cbuff,'(a32,f6.2)') wstr, temp
        call print_log_message(trim(cbuff),5)

        !
        ! Time step
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"reading NIS at position 3 (time step)")
            sz = nis_get_size()
            rbuff(1:1,1:1) = nis_get_real(sz)
            dt = rbuff(1,1)
        else
            call print_error_message(-1,"Missing record at NIS at position 3 (time step)")
        end if
        wstr=adjustl("timeStep    ")
        write(cbuff,'(a32,f9.5)') wstr, dt
        call print_log_message(trim(cbuff),5)

        !
        ! Grid steps
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"reading NIS at position 4 (grid steps)")
            sz = nis_get_size()
            ibuff(1:3,1:1) = nis_get_integer(sz)
            gt = ibuff(1:3,1)
        else
            call print_error_message(-1,"Missing record at NIS at position 4 (grid steps)")
        end if
        wstr = adjustl("gridSteps   ")
        write(cbuff,'(a32,3i4)') wstr, gt
        call print_log_message(trim(cbuff),5)

        !
        ! Grid extent
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"reading NIS at position 5 (grid extent)")
            sz = nis_get_size()
            ibuff(1:3,1:1) = nis_get_integer(sz)
            Nt = ibuff(1:3,1)
            Nte = Nt(1)
        else
            call print_error_message(-1,"Missing record at NIS at position 5 (grid extent)")
        end if
        wstr=adjustl("gridExtent   ")
        write(cbuff,'(a32,3i4)') wstr, Nt
        call print_log_message(trim(cbuff),5)

        !
        ! Calculation of the grid parameters
        grid_Nt = max(gt(1)*Nt(1),gt(2)*Nt(2),gt(3)*Nt(3))
        allocate(grid_t(grid_Nt))                                 ! MARK for deallocation
        do i = 1, grid_Nt
            grid_t(i) = real((i-1),dp)*dt
        end do

        !
        ! save goft ?
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"reading NIS at position 6 (saveGoft)")
            sz = nis_get_size()
            write(cbuff,'(i3)') sz(1)
            write(save_goft,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
        else
            call print_error_message(-1,"Missing record at NIS at position 6 (saveGoft)")
        end if
        ! check allowed values
        call resources_check_values(save_goft,RESOURCES_YES_NO,err)
        wstr=adjustl("saveGoft     ")
        write(cbuff,'(a32,a)') wstr, trim(save_goft)
        call print_log_message(trim(cbuff),5)
        call print_error_message(err,"unknown values for 'saveGoft'")

        !
        ! Save NIS
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"reading NIS at position 7 (save NIS)")
            sz = nis_get_size()
            write(cbuff,'(i3)') sz(1)
            write(save_nis,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
        else
            call print_error_message(-1,"Missing record at NIS at position 7 (save NIS)")
        end if
        ! check allowed values
        call resources_check_values(save_nis,RESOURCES_YES_NO,err)
        wstr=adjustl("saveNIS     ")
        write(cbuff,'(a32,a)') wstr, trim(save_nis)
        call print_log_message(trim(cbuff),5)
        call print_error_message(err,"unknown values for 'saveNIS'")


        !
        ! Parallel
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"reading NIS at position 8 (parallel)")
            sz = nis_get_size()
            write(cbuff,'(i3)') sz(1)
            write(parallel,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
        else
            call print_error_message(-1,"Missing record at NIS at position 8 (parallel)")
        end if
        ! check allowed values
        call resources_check_values(parallel,RESOURCES_YES_NO,err)
        wstr = adjustl("parallel     ")
        write(cbuff,'(a32,a)') wstr, trim(parallel)
        call print_log_message(trim(cbuff),5)
        call print_error_message(err,"unknown values for 'parallel'")

        !
        ! inputFile
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"reading NIS at position 9 (main input file)")
            sz = nis_get_size()
            write(cbuff,'(i3)') sz(1)
            write(inp_file,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
        else
            call print_error_message(-1,"Missing record at NIS at position 9 (main input file)")
        end if
        wstr = adjustl("inputFile     ")
        write(cbuff,'(a32,a)') wstr, trim(inp_file)
        call print_log_message(trim(cbuff),5)

        !
        ! moduleInpFile
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"reading NIS at position 10 (module input file)")
            sz = nis_get_size()
            write(cbuff,'(i3)') sz(1)
            write(mod_inp_file,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
        else
            call print_error_message(-1,"Missing record at NIS at position 10 (module input file)")
        end if
        wstr=adjustl("moduleInpFile     ")
        write(cbuff,'(a32,a)') wstr, trim(mod_inp_file)
        call print_log_message(trim(cbuff),5)

        !
        ! Log Level
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"reading NIS at position 11 (log level)")
            sz = nis_get_size()
            if ((sz(2) == 1).and.(sz(1) == 1)) then
                ibuff(1,1) = nis_get_integer()
                call set_log_level(ibuff(1,1))
            else

            end if
        else
            call print_error_message(-1,"Missing record at NIS at position 11 (log level)")
        end if
        wstr=adjustl("logLevel")
        write(cbuff,'(a32,i2)') wstr, get_log_level()
        call print_log_message(trim(cbuff),5)

        !
        ! outDir
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"reading NIS at position 12 (outputDir)")
            sz = nis_get_size()
            write(cbuff,'(i3)') sz(1)
            write(out_dir,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
        else
            call print_error_message(-1,"Missing record at NIS at position 12 (outputDir)")
        end if
        wstr = adjustl("outputDir")
        write(cbuff,'(a32,a)') wstr, trim(out_dir)
        call print_log_message(trim(cbuff),5)

!        !
!        ! completeResultFile
!        !
!        if (nis_has_next()) then
!            call nis_next(err)
!            call print_error_message(err,"reading NIS at position 12.5 (completeResultFile)")
!            sz = nis_get_size()
!            write(cbuff,'(i3)') sz(1)
!            write(completeResultFile,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
!        else
!            call print_error_message(-1,"Missing record at NIS at position 12.5 (completeResultFile)")
!        end if
!        wstr = adjustl("completeResultFile")
!        write(cbuff,'(a32,a)') wstr, trim(completeResultFile)
!        call print_log_message(trim(cbuff),5)

        !
        ! Number of realizations
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"reading NIS at position 13 (realizations)")
            sz = nis_get_size()
            if ((sz(2) == 1).and.(sz(1) == 1)) then
                N_realizations = nis_get_integer()
            else

            end if
        else
            call print_error_message(-1,"Missing record at NIS at position 13 (realizations)")
        end if
        wstr = adjustl("realizations")
        write(cbuff,'(a32,i5)') wstr, N_realizations
        call print_log_message(trim(cbuff),5)

        !
        ! Frequency of restarts
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"reading NIS at position 14 (restart frequency)")
            sz = nis_get_size()
            if ((sz(2) == 1).and.(sz(1) == 1)) then
                restartFreq = nis_get_integer()
            else

            end if
        else
            call print_error_message(-1,"Missing record at NIS at position 14 (restart frequency)")
        end if
        wstr = adjustl("restartFreq")
        write(cbuff,'(a32,i5)')  wstr,restartFreq
        call print_log_message(trim(cbuff),5)


        !
        ! Outputs
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"reading NIS at position 14.5 (output)")
            sz = nis_get_size()
            write(cbuff,'(i3)') sz(1)
            write(outp,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
            if (outp=="output") then
				wstr = "Reading output specification ..."
            else
	            call nis_next(err)
    	        call print_error_message(err,"output specification expected here")
            end if
        else
            call print_error_message(-1,"Missing record at NIS at position 14.5 (output)")
        end if
        write(cbuff,'(a32,i5)') adjustl(wstr)
        call print_log_message(trim(cbuff),5)
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"reading NIS at position 15 (output nr)")
            sz = nis_get_size()
            if ((sz(2) == 1).and.(sz(1) == 1)) then
                outnr = nis_get_integer()
            else

            end if
        else
            call print_error_message(-1,"Missing record at NIS at position 15 (output nr)")
        end if
        wstr = "number of outputs"
        write(cbuff,'(a32,i5)') adjustl(wstr), outnr
        call print_log_message(trim(cbuff),5)

        do i = 1, outnr

            !
            ! output
            !
            if (nis_has_next()) then
                call nis_next(err)
                call print_error_message(err,"reading NIS at position 16 (output)")
                sz = nis_get_size()
                write(cbuff,'(i3)') sz(1)
                write(outp,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
                call resources_new_output(trim(outp))
            else
                call print_error_message(-1,"Missing record at NIS at position 16 (output)")
            end if
            wstr = "output"
            write(cbuff,'(a32,a)') adjustl(wstr), trim(outp)
            call print_log_message(trim(cbuff),5)

            !
            ! parameters
            !
            if (nis_has_next()) then
                call nis_next(err)
                call print_error_message(err,"reading NIS at position 16.5 (output parameters)")
                sz = nis_get_size()
                write(cbuff,'(i3)') sz(1)
                write(outp,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
                call resources_new_outppar(trim(outp))
            else
                call print_error_message(-1,"Missing record at NIS at position 16.5 (output parameters)")
            end if
            wstr = " -  parameters"
            write(cbuff,'(a32,a)') adjustl(wstr), trim(outp)
            call print_log_message(trim(cbuff),5)

        end do

        !
        ! Tau for tau-dependent projectors
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"reading NIS at position 17 (tau)")
            sz = nis_get_size()
            rbuff(1:1,1:1) = nis_get_real(sz)
            tau_of_projector = rbuff(1,1)
        else
            call print_error_message(-1,"Missing record at NIS at position 17 (tau)")
		end if
        wstr=adjustl("tau    ")
        write(cbuff,'(a32,f6.2)') wstr, tau_of_projector
        call print_log_message(trim(cbuff),5)


        !
        ! RWA frequency
        !
        if (nis_has_next()) then
            call nis_next(err)
            call print_error_message(err,"reading NIS at position 18 (RWA frequency)")
            sz = nis_get_size()
            rbuff(1:1,1:1) = nis_get_real(sz)
            rwa = rbuff(1,1)
        else
            call print_error_message(-1,"Missing record at NIS at position 18 (RWA frequency)")
        end if
        wstr=adjustl("RWA frequency    ")
        write(cbuff,'(a32,f18.12)') wstr, rwa
        call print_log_message(trim(cbuff),5)


        call print_log_message(" ",5)


    end subroutine init_simulation

    !
    ! Returns module name obtained from input
    !
    function get_module_name() result(mname)
        character(len=64) :: mname
        mname = trim(modname)
    end function get_module_name

end module simulate_wizard
