!
! Main body of the NOSE program
!
!
! author: Tomas Mancal
! e-mail: tomas.mancal@mff.cuni.cz
!
! last changed: 2007-01-16
!
program main

    use main_comm

    implicit none

!
! DECLARATION
!

    ! error indicator
    integer            :: err
    integer            :: i,j, bc

    !
    character(len=256) :: char_buff

    real(dp), dimension(:,:), pointer :: ah
!
! BODY
!

    !
    ! Initializing parallel processing
    !
    call init_parallel()


    call init_allocation_counting()

    !
    ! Initialization of timing and logging services
    !


    call init_std_io(parallel_id)
    call init_timing()

    call init_resources()


    !
    ! NOSE initialization
    !
    call print_header_parallel()
    call warning_parallel()


    !
    ! Simple message among the processes
    !
    call parallel_communication_test()


    !
    ! Reading input (master reads and distributes)
    !
    call start_section_parallel("Reading the NIS")
    !
    call init_nis()
    call nis_set_mode(NIS_MODE_STD)            ! sets mode to standars i/o
    call nis_read()                            ! read from the standard (std) i/o
    call nis_distribute()                      ! send info to others

    call print_log_message("NIS received",1)   ! respond to tcl script
    call flush()

    !read(*,*) char_buff                        ! wait for answer

    !
    call next_section_parallel("sm")


    !
    call warning_output()


    !
    ! Initialization of the simulation
    !
    call start_section_parallel("Initializing simulation")
    !
    call init_simulation(err)                   ! reading simulation data from NIS
    call print_error_message(err,"reading simulation input")
    !
    call next_section_parallel("sm")
    !
    ! End of initialization
    !


    !
    ! Reading module input
    !
    call start_section_parallel("Reading input")
    !
    ! obtain the module name
    !modname = get_module_name()
    call read_selected_input(err) !trim(modname),err)                 ! reading module input from NIS
    write(char_buff,'(a,a)') "reading module input in: ", trim(modname)
    call print_error_message(err,trim(char_buff))
    call resources_after_input()
    !
    call next_section_parallel("sm")
    !
    ! End of input reading
    !

    !
    ! Save the NIS from all processors to make mutual check
    !
    if (trim(save_nis) == "yes") then
        write(char_buff,'(a,i1,a)') "nis_logfile_",parallel_id,".log"
        call nis_open_file(50+parallel_id,trim(char_buff))
        call nis_set_mode(NIS_MODE_FILE)
        call nis_send()
        call nis_close_file()
    end if



    ! these things will be set somewhere else
    parallel_nr_uspr       = 35
    parallel_nr_runs_total = N_realizations

    !
    ! Distribute seeds for random generators
    !
    call nis_distribute_seeds()


    !
    ! Everything that has to be done before calculating realizations
    !
    call prepare_for_all()


    !***************
    !  Start cycle
    !***************
    call start_section_parallel("Simulating")
    N_since_saved = 0
    do i_realization = 1, N_realizations_local

        N_since_saved = N_since_saved + 1

        !
        !
        ! Calculation of quantities needed for simulation
        !
        call prepare_for_one(err)
        call print_error_message(err,"preparing for simulation")
        !
        ! End
        !


        !
        ! Simulation
        !
        call do_module_work(trim(modname),err)
        if (err < 0) then
            print *, "Error while simulating"
            stop
        end if
        !
        ! End of simulation
        !


        if ((N_since_saved*parallel_nr >= restartFreq).or.(i_realization == N_realizations_local)) then

            write(char_buff,'(a,i5,a)') "Restart point; total of ", &
              i_realization*parallel_nr, " realizations finished"
            call print_log_message(trim(char_buff),5)
            bc = get_byte_count()
            write(char_buff,'(a,f10.4,a)') "currently using ", &
             real(bc)/real(1024**2), " MB of memory"
            call print_log_message(trim(char_buff),5)

            call save_state_for_restart()

            N_since_saved = 0

            call flush()

        end if



    end do
    call next_section_parallel("smh")

!    call parallel_synchronize()


    !
    ! Collect data (send to the main process of wait for them
    ! if you are a main process)
    !
!    call start_section_parallel("Waiting for other processes to deliver data")
!    call collect_module_data(trim(modname),err)
!    if (err < 0) then
!        print *, "Error collecting data"
!        stop
!    end if
    !
    ! End of collecting data
    !



    bc = get_byte_count()
    write(char_buff,'(a,f10.4,a)') "currently using ", &
        real(bc)/real(1024**2), " MB of memory"
    call print_log_message(trim(char_buff),5)

    !********************************
    ! End cycle
    !********************************

    !
    ! Cleaning after work
    !
    call clean_module(err)
    if (err < 0) then
        print *, "Error cleaning memory"
        stop
    end if

    call clean_resources()
    !
    call next_section_parallel("smh")
    !
    ! End of cleaning
    !

    !
    ! Finish the parallel work
    !
    call start_section_parallel("Closing MPI")
    !
    call finalize_parallel(err)
    if (err < 0) then
        print *, "Error in unspliting"
        stop
    end if
    !
    call next_section_parallel("sm")
    !
    ! End
    !

    call clean_nis()
    !
    call finish_sections_parallel("smhd")
    call print_final_parallel()



!
! END BODY
!
contains

    subroutine save_state_for_restart

    	call parallel_synchronize()
		call start_section_parallel("Waiting for other processes to deliver data")
    	call collect_module_data(trim(modname),err)
    	if (err < 0) then
        	print *, "Error collecting data"
        	stop
    	end if

        call print_log_message("... state saved",5)

    end subroutine save_state_for_restart

end program main
!
! End of the program
!
