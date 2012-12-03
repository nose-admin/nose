!
! Input Wizard
!
! Reads input of all modules registred in the NOSE program
!
! Author:       Tomas Mancal
! E-mail:       tomas.mancal@mff.cuni.cz
! Last changed: 2007-01-11
!
module input_wizard

    use simulate_wizard
    use input_tdpt3
    use input_qme
    use input_montecarlo

    !use data_list

    implicit none

    ! Error information
    character(len=32), parameter, private :: file_name = "input_wizard.f90"


contains
    !
    ! Reads input for a given module
    !
    subroutine read_selected_input(err) !modname,err)
        !character(len=*), intent(in) :: modname
        integer, intent(out)          :: err

        err = 0

        call read_input(err)

        !
        ! Module TDPT-3
        !

        if (trim(modname) == "TDPT-3") then

            !
            ! reads special module input
            !
            call read_input_tdpt3(err)


        else if (trim(modname) == "QME") then

            !
            ! reads special module input
            !
            call read_input_qme(err)


		else if (trim(modname) == "MC") then

            !
            ! reads special module input
            !
            call read_input_montecarlo(err)

		else if (trim(modname) == "TEST") then

            !
            ! reads special module input
            !
            !call read_input_test(err)

        else if (trim(modname) == "MOLC") then

        else

            print *, "Unknown module name ", modname
            err = -1

        end if

        !
        ! End of module TDPT-3
        !

    end subroutine read_selected_input

    !
    ! Reads input common to all modules in NOSE
    !
    subroutine read_input(err)
        integer, intent(out)              :: err
        integer, dimension(2)             :: sz
        real(dp), dimension(:), pointer   :: rbuff_v
        integer                              :: ibuff
        integer, dimension(:), pointer    :: ibuff_v
        real(dp), dimension(:,:), pointer :: rbuff_m
        real(dp), dimension(:,:,:), pointer :: rbuff_m2
        real(dp)                          :: rbuff_s
        character(len=512) :: cbuff
        character(len=256) :: caux
        integer :: i


        !type(list_dp) :: param_list
        !type(list_element_dp), pointer :: last_el

        err = 0

        call init_resources()

        !************************************************************************
        ! Reading the input file data
        !************************************************************************

        do
            !
            ! reads one of : SECTION_BLOCKS, SECTION_INTER_BLOCK SECTION_GOFTS
            !
            if (nis_has_next()) then
                call nis_next(err)
                call print_error_message(err,"reading NIS at the section separator")
                sz = nis_get_size()
                write(cbuff,'(i3)') sz(1)
                write(caux,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
            else
                call print_error_message(-1,"Missing record at NIS (section separator")
            end if
            write(cbuff,'(a,a)') "section separator -  ", trim(caux)
            call print_log_message(trim(cbuff),7)


            if (trim(caux) == "SECTION_BLOCKS") then

                do   ! block loop

                    !
                    ! reads the BLOCK_CONTINUE or the BLOCK_END string
                    !
                    if (nis_has_next()) then
                        call nis_next(err)
                        call print_error_message(err,"reading NIS at block separator")
                        sz = nis_get_size()
                        write(cbuff,'(i3)') sz(1)
                        write(caux,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
                    else
                        call print_error_message(-1,"Missing record at NIS (block separator")
                    end if
                    write(cbuff,'(a,a)') "block separator -  ", trim(caux)
                    call print_log_message(trim(cbuff),7)

                    if (trim(caux) == "BLOCK_END") exit

                    if (trim(caux) == "BLOCK_CONTINUE") then

                        call resources_new_block()

                        !
                        ! Reading block id
                        !
                        if (nis_has_next()) then
                            call nis_next(err)
                            call print_error_message(err,"reading NIS at block id")
                            sz = nis_get_size()
                            ibuff = nis_get_integer()
                            current_s_block%id = ibuff
                            !current_e_block%id = ibuff
                        else
                            call print_error_message(-1,"Missing record at NIS (block id")
                        end if
                        write(cbuff,'(a,a)') "block id -  ", ibuff

                        !
                        ! Site Energies
                        !
                        if (nis_has_next()) then
                            call nis_next(err)
                            call print_error_message(err,"reading site energies from NIS")
                            sz = nis_get_size()

                            ! number of one exciton states
                            allocate(current_s_block%N1)
                            !allocate(current_e_block%N1)
                            current_s_block%N1 = sz(1)
                            !current_e_block%N1 = sz(1)

                            allocate(rbuff_v(sz(1)))
                            rbuff_v = nis_get_real(sz(1))
                            current_s_block%en => rbuff_v(1:sz(1))
                            allocate(current_s_block%en_orig(1:sz(1)))
                            current_s_block%en_orig = current_s_block%en


                        else
                            call print_error_message(-1,"Missing record at NIS (site energies)")
                        end if
                        write(cbuff,'(a,3i7)') "site energies read succesfully "
                        call print_log_message(trim(cbuff),7)

                        !
                        ! Positions
                        !
                        if (nis_has_next()) then
                            call nis_next(err)
                            call print_error_message(err,"reading dipole coordinates from NIS")
                            sz = nis_get_size()
                            if (sz(1) /= current_s_block%N1) then
                                call print_error_message(-1,"dipole coordinates from NIS - wrong number of dipoles")
                            end if
                            allocate(rbuff_m(sz(1),sz(2)))
                            rbuff_m = nis_get_real(sz)
                            current_s_block%rr => rbuff_m
                            !if (parallel_id == 0) then
                            !    do i = 1, current_s_block%N1
                            !        print *, current_s_block%rr(i,:)
                            !    end do
                            !end if
                        else
                            call print_error_message(-1,"Missing record at NIS (dipole coordinates)")
                        end if


                        !
                        ! Dipole monents
                        !

                        ! x component
                        if (nis_has_next()) then
                            call nis_next(err)
                            call print_error_message(err,"reading dipole x-component from NIS")
                            sz = nis_get_size()
                            if (sz(1) /= current_s_block%N1) then
                                call print_error_message(-1,"dipole x-component from NIS - wrong number of dipoles")
                            end if
                            allocate(rbuff_v(sz(1)))
                            rbuff_v = nis_get_real(sz(1))
                            current_s_block%dx => rbuff_v(1:sz(1))
                        else
                            call print_error_message(-1,"Missing record at NIS (dipole x-component)")
                        end if

                        ! y component
                        if (nis_has_next()) then
                            call nis_next(err)
                            call print_error_message(err,"reading dipole y-component from NIS")
                            sz = nis_get_size()
                            if (sz(1) /= current_s_block%N1) then
                                call print_error_message(-1,"dipole y-component from NIS - wrong number of dipoles")
                            end if
                            allocate(rbuff_v(sz(1)))
                            rbuff_v = nis_get_real(sz(1))
                            current_s_block%dy => rbuff_v(1:sz(1))
                        else
                            call print_error_message(-1,"Missing record at NIS (dipole x-component)")
                        end if

                        ! z component
                        if (nis_has_next()) then
                            call nis_next(err)
                            call print_error_message(err,"reading dipole z-component from NIS")
                            sz = nis_get_size()
                            if (sz(1) /= current_s_block%N1) then
                                call print_error_message(-1,"dipole z-component from NIS - wrong number of dipoles")
                            end if
                            allocate(rbuff_v(sz(1)))
                            rbuff_v = nis_get_real(sz(1))
                            current_s_block%dz => rbuff_v(1:sz(1))
                        else
                            call print_error_message(-1,"Missing record at NIS (dipole x-component)")
                        end if
                        write(cbuff,'(a,3i7)') "transition dipole moments read succesfully "
                        call print_log_message(trim(cbuff),7)


                        ! normalization of dipoles
                        do i = 1, current_s_block%N1
                            rbuff_s = sqrt(current_s_block%dx(i)**2 + current_s_block%dy(i)**2 &
                            + current_s_block%dz(i)**2)
                            !rbuff_s = 1.0_dp
                            current_s_block%dx(i) = current_s_block%dx(i)/rbuff_s
                            current_s_block%dy(i) = current_s_block%dy(i)/rbuff_s
                            current_s_block%dz(i) = current_s_block%dz(i)/rbuff_s
                        end do


                        !
                        ! Dipole length
                        !
                        if (nis_has_next()) then
                            call nis_next(err)
                            call print_error_message(err,"reading dipole lengths from NIS")
                            sz = nis_get_size()
                            if (sz(1) /= current_s_block%N1) then
                                call print_error_message(-1,"dipole length from NIS - wrong number of dipoles")
                            end if
                            allocate(rbuff_v(sz(1)))
                            rbuff_v = nis_get_real(sz(1))
                            !print *, rbuff_v
                            current_s_block%dx(1:sz(1)) = current_s_block%dx(1:sz(1))*rbuff_v(1:sz(1))
                            current_s_block%dy(1:sz(1)) = current_s_block%dy(1:sz(1))*rbuff_v(1:sz(1))
                            current_s_block%dz(1:sz(1)) = current_s_block%dz(1:sz(1))*rbuff_v(1:sz(1))

                        else
                            call print_error_message(-1,"Missing record at NIS (dipole length)")
                        end if
                        write(cbuff,'(a,3i7)') "transition dipole moment lengths read succesfully "
                        call print_log_message(trim(cbuff),7)

                        !
                        ! Coupling
                        !
                        if (nis_has_next()) then
                            call nis_next(err)
                            call print_error_message(err,"reading coupling energies from NIS")
                            sz = nis_get_size()
                            if ((sz(1) /= current_s_block%N1).or.(sz(2) /= current_s_block%N1)) then
                                write(cbuff,'(a,2i3,a,i3)') "coupling energies - wrong rank of the matrix; is ", &
                                       sz(1), sz(2), " should be ", current_s_block%N1
                                call print_error_message(-1,trim(cbuff))
                            end if
                            allocate(rbuff_m(sz(1),sz(2)))
                            rbuff_m = nis_get_real(sz)
                            current_s_block%J => rbuff_m
                        else
                            call print_error_message(-1,"Missing record at NIS (coupling energies)")
                        end if
                        write(cbuff,'(a)') "coupling energies read succesfully "
                        call print_log_message(trim(cbuff),7)

                        ! gofts
                        if (nis_has_next()) then
                            call nis_next(err)
                            call print_error_message(err,"reading goft index from NIS")
                            sz = nis_get_size()
                            if (sz(1) /= current_s_block%N1) then
                                call print_error_message(-1,"goft index from NIS - wrong size")
                            end if
                            allocate(ibuff_v(sz(1)))
                            ibuff_v = nis_get_integer(sz(1))
                            current_s_block%gindex => ibuff_v(1:sz(1))
                        else
                            call print_error_message(-1,"Missing record at NIS (goft index)")
                        end if
                        write(cbuff,'(a,3i7)') "gindex read succesfully "
                        call print_log_message(trim(cbuff),7)

                        ! disorder indices
                        if (nis_has_next()) then
                            call nis_next(err)
                            call print_error_message(err,"reading disorder index from NIS")
                            sz = nis_get_size()
                            if (sz(1) /= current_s_block%N1) then
                                call print_error_message(-1,"disorder index from NIS - wrong size")
                            end if
                            allocate(ibuff_v(sz(1)))
                            ibuff_v = nis_get_integer(sz(1))
                            current_s_block%dindex => ibuff_v(1:sz(1))
                        else
                            call print_error_message(-1,"Missing record at NIS (disorder index)")
                        end if
                        write(cbuff,'(a,3i7)') "dindex read succesfully "
                        call print_log_message(trim(cbuff),7)

                        ! disorder width
                        if (nis_has_next()) then
                            call nis_next(err)
                            call print_error_message(err,"reading disorder widths from NIS")
                            sz = nis_get_size()
                            if ((sz(1) /= current_s_block%N1).or.(sz(2) /= current_s_block%N1)) then
                                write(cbuff,'(a,2i3,a,i3)') "disorder widths - wrong rank of the matrix; is ", &
                                       sz(1), sz(2), " should be ", current_s_block%N1
                                call print_error_message(-1,trim(cbuff))
                            end if
                            allocate(rbuff_m(sz(1),sz(2)))
                            rbuff_m = nis_get_real(sz)
                            current_s_block%dwidth => rbuff_m
                        else
                            call print_error_message(-1,"Missing record at NIS (disorder widths)")
                        end if
                        write(cbuff,'(a)') "disorder width read succesfully "
                        call print_log_message(trim(cbuff),7)

                       ! rotational strengths
                        if (nis_has_next()) then
                            call nis_next(err)
                            call print_error_message(err,"reading rotational strengths from NIS")
                            sz = nis_get_size()
                            if (sz(1) /= current_s_block%N1) then
                                call print_error_message(-1,"rotational strengths from NIS - wrong number of dipoles")
                            end if
                            allocate(rbuff_v(sz(1)))
                            rbuff_v = nis_get_real(sz(1))
                            current_s_block%rs => rbuff_v(1:sz(1))
                        else
                            call print_error_message(-1,"Missing record at NIS (rotational strengths)")
                        end if


                    else

                        call print_error_message(-1,"wrong block separator")

                    end if

                end do ! end block loop

            else if (trim(caux) == "SECTION_PERIODIC") then


                    !
                    ! Reading presence of PERIODIC section
                    !
                    if (nis_has_next()) then
                     	call nis_next(err)
                        call print_error_message(err,"reading NIS at periodic")
                        sz = nis_get_size()
                        ibuff = nis_get_integer()
                        current_s_block%periodic = ibuff

                    else
                        call print_error_message(-1,"Missing record at NIS (periodic)")
                    end if
                    write(cbuff,'(a,a)') "periodic -  ", ibuff

                    if (current_s_block%periodic == 1) then

                    	!
                    	! Reading periodicity dimension
                    	!
                    	if (nis_has_next()) then
                     		call nis_next(err)
                        	call print_error_message(err,"reading NIS at periodicity dimension")
                        	sz = nis_get_size()
                        	ibuff = nis_get_integer()
                        	current_s_block%dimension = ibuff

                    	else
                        	call print_error_message(-1,"Missing record at NIS (periodicity dimension)")
                    	end if
                    	write(cbuff,'(a,a)') "periodicity dimension -  ", ibuff

                    	!
                    	! Reading periodicity
                    	!
                    	if (nis_has_next()) then
                     		call nis_next(err)
                        	call print_error_message(err,"reading NIS at periodicity")
                        	sz = nis_get_size()
                        	ibuff = nis_get_integer()
                        	current_s_block%periodicity = ibuff

                    	else
                        	call print_error_message(-1,"Missing record at NIS (periodicity)")
                    	end if
                    	write(cbuff,'(a,a)') "periodicity -  ", ibuff


                        ! R1 vector
                        if (nis_has_next()) then
                            call nis_next(err)
                            call print_error_message(err,"reading R1 from NIS")
                            sz = nis_get_size()
                            if (sz(1) /= 3) then
                                call print_error_message(-1,"R1 from NIS - wrong number of component")
                            end if
                            allocate(rbuff_v(sz(1)))
                            rbuff_v = nis_get_real(sz(1))
                            current_s_block%cell_vec1 => rbuff_v(1:sz(1))
                        else
                            call print_error_message(-1,"Missing record at NIS (R1)")
                        end if

                        ! R2 vector
                        if (nis_has_next()) then
                            call nis_next(err)
                            call print_error_message(err,"reading R2 from NIS")
                            sz = nis_get_size()
                            if (sz(1) /= 3) then
                                call print_error_message(-1,"R2 from NIS - wrong number of component")
                            end if
                            allocate(rbuff_v(sz(1)))
                            rbuff_v = nis_get_real(sz(1))
                            current_s_block%cell_vec2 => rbuff_v(1:sz(1))
                        else
                            call print_error_message(-1,"Missing record at NIS (R2)")
                        end if

                        ! R3 vector
                        if (nis_has_next()) then
                            call nis_next(err)
                            call print_error_message(err,"reading R3 from NIS")
                            sz = nis_get_size()
                            if (sz(1) /= 3) then
                                call print_error_message(-1,"R3 from NIS - wrong number of component")
                            end if
                            allocate(rbuff_v(sz(1)))
                            rbuff_v = nis_get_real(sz(1))
                            current_s_block%cell_vec3 => rbuff_v(1:sz(1))
                        else
                            call print_error_message(-1,"Missing record at NIS (R3)")
                        end if

						if (current_s_block%dimension == 1) then
	                        !
    	                    ! Neighbor cell coupling
        	                !
            	            if (nis_has_next()) then
                	            call nis_next(err)
                    	        call print_error_message(err,"reading neighbor coupling energies from NIS")
                        	    sz = nis_get_size()
                            	if ((sz(1) /= current_s_block%N1).or.(sz(2) /= current_s_block%N1)) then
                                	write(cbuff,'(a,2i3,a,i3)') "neighbor coupling energies - wrong rank of the matrix; is ", &
                                       sz(1), sz(2), " should be ", current_s_block%N1
                                	call print_error_message(-1,trim(cbuff))
                            	end if
                            	allocate(rbuff_m2(1,sz(1),sz(2)))
                            	rbuff_m2(1,:,:) = nis_get_real(sz)
                            	current_s_block%Jcell => rbuff_m2
                        	else
                            	call print_error_message(-1,"Missing record at NIS (coupling energies)")
                        	end if
                        	write(cbuff,'(a)') "neighbor coupling energies read succesfully "
                        	call print_log_message(trim(cbuff),7)

						else if (current_s_block%dimension == 2) then
	                        !
    	                    ! Neighbor cell coupling (2D)
        	                !
            	            if (nis_has_next()) then
                	            call nis_next(err)
                    	        call print_error_message(err,"reading neighbor (2D) coupling energies from NIS")
                        	    sz = nis_get_size()
                            	if ((sz(1) /= current_s_block%N1).or.(sz(2) /= current_s_block%N1)) then
                                	write(cbuff,'(a,2i3,a,i3)') "neighbor coupling (2D) energies - wrong rank of the matrix; is ", &
                                       sz(1), sz(2), " should be ", current_s_block%N1
                                	call print_error_message(-1,trim(cbuff))
                            	end if
                            	allocate(rbuff_m2(4,sz(1),sz(2)))
                            	rbuff_m2(1,:,:) = nis_get_real(sz)
                        	else
                            	call print_error_message(-1,"Missing record at NIS (neighbor (2D) coupling energies)")
                        	end if
                        	write(cbuff,'(a)') "neighbor coupling energies read succesfully "
                        	call print_log_message(trim(cbuff),7)

            	            if (nis_has_next()) then
                	            call nis_next(err)
                    	        call print_error_message(err,"reading neighbor (2D) coupling energies from NIS")
                        	    sz = nis_get_size()
                            	if ((sz(1) /= current_s_block%N1).or.(sz(2) /= current_s_block%N1)) then
                                	write(cbuff,'(a,2i3,a,i3)') "neighbor coupling (2D) energies - wrong rank of the matrix; is ", &
                                       sz(1), sz(2), " should be ", current_s_block%N1
                                	call print_error_message(-1,trim(cbuff))
                            	end if
                            	rbuff_m2(2,:,:) = nis_get_real(sz)
                        	else
                            	call print_error_message(-1,"Missing record at NIS (neighbor (2D) coupling energies)")
                        	end if
                        	write(cbuff,'(a)') "neighbor coupling energies read succesfully "
                        	call print_log_message(trim(cbuff),7)

            	            if (nis_has_next()) then
                	            call nis_next(err)
                    	        call print_error_message(err,"reading neighbor (2D) coupling energies from NIS")
                        	    sz = nis_get_size()
                            	if ((sz(1) /= current_s_block%N1).or.(sz(2) /= current_s_block%N1)) then
                                	write(cbuff,'(a,2i3,a,i3)') "neighbor coupling (2D) energies - wrong rank of the matrix; is ", &
                                       sz(1), sz(2), " should be ", current_s_block%N1
                                	call print_error_message(-1,trim(cbuff))
                            	end if
                            	rbuff_m2(3,:,:) = nis_get_real(sz)
                        	else
                            	call print_error_message(-1,"Missing record at NIS (neighbor (2D) coupling energies)")
                        	end if
                        	write(cbuff,'(a)') "neighbor coupling energies read succesfully "
                        	call print_log_message(trim(cbuff),7)

            	            if (nis_has_next()) then
                	            call nis_next(err)
                    	        call print_error_message(err,"reading neighbor (2D) coupling energies from NIS")
                        	    sz = nis_get_size()
                            	if ((sz(1) /= current_s_block%N1).or.(sz(2) /= current_s_block%N1)) then
                                	write(cbuff,'(a,2i3,a,i3)') "neighbor coupling (2D) energies - wrong rank of the matrix; is ", &
                                       sz(1), sz(2), " should be ", current_s_block%N1
                                	call print_error_message(-1,trim(cbuff))
                            	end if
                            	rbuff_m2(4,:,:) = nis_get_real(sz)
                        	else
                            	call print_error_message(-1,"Missing record at NIS (neighbor (2D) coupling energies)")
                        	end if
                        	write(cbuff,'(a)') "neighbor coupling energies read succesfully "
                        	call print_log_message(trim(cbuff),7)

                            current_s_block%Jcell => rbuff_m2

						else

							call print_error_message(-1,"Dimensions other then 1 and 2 are not yet supported")

						end if

                    end if

                    !
                    ! reads untill the PERIODIC_END string
                    !
                    if (nis_has_next()) then
                        call nis_next(err)
                        call print_error_message(err,"reading NIS at PERODIC")
                        sz = nis_get_size()
                        write(cbuff,'(i3)') sz(1)
                        write(caux,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
                    else
                        call print_error_message(-1,"Missing record at NIS (PERIODIC)")
                    end if
                    write(cbuff,'(a,a)') "PERIODIC  ", trim(caux)
                    call print_log_message(trim(cbuff),7)

                    if (trim(caux) /= "PERIODIC_END") then
                    	call print_error_message(-1,"End of PERIODIC section expected");
                    end if


            else if (trim(caux) == "SECTION_INTER_BLOCKS") then

                 do   ! inter block loop

                    !
                    ! reads the INTER_BLOCK_CONTINUE or the INTER_BLOCK_END string
                    !
                    if (nis_has_next()) then
                        call nis_next(err)
                        call print_error_message(err,"reading NIS at goft separator")
                        sz = nis_get_size()
                        write(cbuff,'(i3)') sz(1)
                        write(caux,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
                    else
                        call print_error_message(-1,"Missing record at NIS (goft separator")
                    end if
                    write(cbuff,'(a,a)') "goft separator -  ", trim(caux)
                    call print_log_message(trim(cbuff),7)

                    if (trim(caux) == "INTER_BLOCK_END") exit

                    if (trim(caux) == "INTER_BLOCK_CONTINUE") then

                        call resources_new_inter_block()

                        !
                        ! Reading inter block id1
                        !
                        if (nis_has_next()) then
                            call nis_next(err)
                            call print_error_message(err,"reading NIS at inter_block id1")
                            sz = nis_get_size()
                            ibuff = nis_get_integer()

                            current_i_block%id = ibuff

                        else
                            call print_error_message(-1,"Missing record at NIS (inter_block id1")
                        end if
                        write(cbuff,'(a,i3)') " id1 -  ", ibuff
                        call print_log_message(trim(cbuff),7)

                        !
                        ! Reading inter block id2
                        !
                        if (nis_has_next()) then
                            call nis_next(err)
                            call print_error_message(err,"reading NIS at inter_block id2")
                            sz = nis_get_size()
                            ibuff = nis_get_integer()

                            !current_i_block%id = ibuff

                        else
                            call print_error_message(-1,"Missing record at NIS (inter_block id2")
                        end if
                        write(cbuff,'(a,i3)') " id2 -  ", ibuff
                        call print_log_message(trim(cbuff),7)


                        !
                        ! Inter Coupling
                        !
                        if (nis_has_next()) then
                            call nis_next(err)
                            call print_error_message(err,"reading inter block coupling energies from NIS")
                            sz = nis_get_size()
                            !if ((sz(1) /= current_s_block%N1).or.(sz(2) /= current_s_block%N1)) then
                            !    write(cbuff,'(a,2i3,a,i3)') "coupling energies - wrong rank of the matrix; is ", &
                            !           sz(1), sz(2), " should be ", current_s_block%N1
                            !    call print_error_message(-1,trim(cbuff))
                            !end if
                            allocate(rbuff_m(sz(1),sz(2)))
                            rbuff_m = nis_get_real(sz)
                            current_i_block%J => rbuff_m

                            !print *, current_i_block%J
                        else
                            call print_error_message(-1,"Missing record at NIS (inter coupling energies)")
                        end if
                        write(cbuff,'(a)') "inter coupling energies read succesfully "
                        call print_log_message(trim(cbuff),7)



                    end if

                 end do

            else if (trim(caux) == "SECTION_GOFTS") then

                 do   ! goft loop

                    !
                    ! reads the GOFT_CONTINUE or the GOFT_END string
                    !
                    if (nis_has_next()) then
                        call nis_next(err)
                        call print_error_message(err,"reading NIS at goft separator")
                        sz = nis_get_size()
                        write(cbuff,'(i3)') sz(1)
                        write(caux,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
                    else
                        call print_error_message(-1,"Missing record at NIS (goft separator")
                    end if
                    write(cbuff,'(a,a)') "goft separator -  ", trim(caux)
                    call print_log_message(trim(cbuff),7)

                    if (trim(caux) == "GOFT_END") exit

                    if (trim(caux) == "GOFT_CONTINUE") then

                        call resources_new_goft()

                        !
                        ! Reading goft id
                        !
                        if (nis_has_next()) then
                            call nis_next(err)
                            call print_error_message(err,"reading NIS at goft id")
                            sz = nis_get_size()
                            ibuff = nis_get_integer()

                            current_s_goft%id = ibuff

                        else
                            call print_error_message(-1,"Missing record at NIS (goft id")
                        end if
                        write(cbuff,'(a,i3)') "goft id -  ", ibuff
                        call print_log_message(trim(cbuff),7)

                        do

                            !
                            ! reads the g(t) mode type the MODE_END string
                            !
                            ! types:
                            !       BROWNIAN
                            !       BROWNIAN_GENERAL
                            !       BROWNIAN_UNDERDAMPED
                            !       DELTA
                            !       GAUSSIAN
                            !
                            if (nis_has_next()) then
                                call nis_next(err)
                                call print_error_message(err,"reading NIS at goft type")
                                sz = nis_get_size()
                                write(cbuff,'(i3)') sz(1)
                                write(caux,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
                            else
                                call print_error_message(-1,"Missing record at NIS (goft type")
                            end if
                            write(cbuff,'(a,a)') "goft type -  ", trim(caux)
                            call print_log_message(trim(cbuff),7)

                            if (trim(caux) == "MODE_END") exit

                            if (trim(caux) == "BROWNIAN") then

                                current_s_goft%nr_modes = current_s_goft%nr_modes + 1
                                current_s_goft%types(current_s_goft%nr_modes) = "BROWNIAN"

                                !
                                ! Brownian parameters
                                !
                                if (nis_has_next()) then
                                    call nis_next(err)
                                    call print_error_message(err,"reading Brownian mode from NIS")
                                    sz = nis_get_size()
                                    if ((sz(1) /= 2).or.(sz(2) /= 1)) then
                                        call print_error_message(-1,"Brownian mode from NIS - wrong rank")
                                    end if
                                    allocate(rbuff_v(2))
                                    rbuff_v = nis_get_real(2)
                                    current_s_goft%params(1:2,current_s_goft%nr_modes) = rbuff_v
                                else
                                    call print_error_message(-1,"Missing record at NIS (Brownian mode)")
                                end if
                                write(cbuff,'(a,2f18.12)') "Brownian mode parameters ", rbuff_v
                                call print_log_message(trim(cbuff),7)

                            else if (trim(caux) == "GAUSSIAN") then

                                current_s_goft%nr_modes = current_s_goft%nr_modes + 1
                                current_s_goft%types(current_s_goft%nr_modes) = "GAUSSIAN"


                                !
                                ! Brownian parameters
                                !
                                if (nis_has_next()) then
                                    call nis_next(err)
                                    call print_error_message(err,"reading Gaussian mode from NIS")
                                    sz = nis_get_size()
                                    if ((sz(1) /= 2).or.(sz(2) /= 1)) then
                                        call print_error_message(-1,"Gaussian mode from NIS - wrong rank")
                                    end if
                                    allocate(rbuff_v(2))
                                    rbuff_v = nis_get_real(2)
                                    current_s_goft%params(1:2,current_s_goft%nr_modes) = rbuff_v
                                else
                                    call print_error_message(-1,"Missing record at NIS (Gaussian mode)")
                                end if
                                write(cbuff,'(a,2f18.12)') "Gaussian mode parameters ", rbuff_v
                                call print_log_message(trim(cbuff),7)

                            else if (trim(caux) == "DELTA") then

                                current_s_goft%nr_modes = current_s_goft%nr_modes + 1
                                current_s_goft%types(current_s_goft%nr_modes) = "DELTA"

                                !
                                ! Delta mode parameters
                                !
                                if (nis_has_next()) then
                                    call nis_next(err)
                                    call print_error_message(err,"reading Delta mode from NIS")
                                    sz = nis_get_size()
                                    if ((sz(1) /= 1).or.(sz(2) /= 1)) then
                                        call print_error_message(-1,"Delta mode from NIS - wrong rank")
                                    end if
                                    allocate(rbuff_v(1))
                                    rbuff_v = nis_get_real(1)
                                    current_s_goft%params(1:1,current_s_goft%nr_modes) = rbuff_v
                                else
                                    call print_error_message(-1,"Missing record at NIS (Brownian mode)")
                                end if
                                write(cbuff,'(a,2f18.12)') "Delta mode parameters ", rbuff_v
                                call print_log_message(trim(cbuff),7)

                            else if (trim(caux) == "BROWNIAN_GENERAL") then

                                current_s_goft%nr_modes = current_s_goft%nr_modes + 1
                                current_s_goft%types(current_s_goft%nr_modes) = "BROWNIAN_GENERAL"

                                !
                                ! Brownian-general parameters
                                !
                                if (nis_has_next()) then
                                    call nis_next(err)
                                    call print_error_message(err,"reading Brownian-general mode from NIS")
                                    sz = nis_get_size()
                                    if ((sz(1) /= 3).or.(sz(2) /= 1)) then
                                        call print_error_message(-1,"Brownian-general mode from NIS - wrong rank")
                                    end if
                                    allocate(rbuff_v(3))
                                    rbuff_v = nis_get_real(3)
                                    current_s_goft%params(1:3,current_s_goft%nr_modes) = rbuff_v
                                else
                                    call print_error_message(-1,"Missing record at NIS (Brownian-general mode)")
                                end if
                                write(cbuff,'(a,3f18.12)') "Brownian-general mode parameters ", rbuff_v
                                call print_log_message(trim(cbuff),7)

                            else if (trim(caux) == "BROWNIAN_UNDERDAMPED") then

                                current_s_goft%nr_modes = current_s_goft%nr_modes + 1
                                current_s_goft%types(current_s_goft%nr_modes) = "BROWNIAN_UNDERDAMPED"

                                !
                                ! Brownian-underdamped parameters
                                !
                                if (nis_has_next()) then
                                    call nis_next(err)
                                    call print_error_message(err,"reading Brownian-underdamped mode from NIS")
                                    sz = nis_get_size()
                                    if ((sz(1) /= 2).or.(sz(2) /= 1)) then
                                        call print_error_message(-1,"Brownian-underdamped mode from NIS - wrong rank")
                                    end if
                                    allocate(rbuff_v(2))
                                    rbuff_v = nis_get_real(2)
                                    current_s_goft%params(1:2,current_s_goft%nr_modes) = rbuff_v
                                else
                                    call print_error_message(-1,"Missing record at NIS (Brownian-general mode)")
                                end if
                                write(cbuff,'(a,2f18.12)') "Brownian-underdamped mode parameters ", rbuff_v
                                call print_log_message(trim(cbuff),7)

                            else

                                call print_error_message(-1,"unexpected g(t) type")

                            end if

                        end do

                    end if

                end do ! END goft loop

                !call list_show(param_list)
                !call list_destroy(param_list)

            else if (trim(caux) == "SECTION_INTER_BLOCK") then


            else if (trim(caux) == "SECTION_DM_INITIAL_CONDITION") then

            	allocate(ini_dm(current_s_block%N1,current_s_block%N1))

            	ini_dm = 0.0d0


                    !
                    ! Basis type
                    !
                    if (nis_has_next()) then
                    	call nis_next(err)
                        call print_error_message(err,"reading NIS at basis type")
                        sz = nis_get_size()
                        ibuff = nis_get_integer()
                        ini_basis_type = ibuff

                    else
                    	call print_error_message(-1,"Missing record at NIS (basis type)")
                    end if
                    write(cbuff,'(a,i5)') "basis type -  ", ibuff

                    !
                    ! Number of states
                    !
                    if (nis_has_next()) then
                    	call nis_next(err)
                        call print_error_message(err,"reading NIS at number of states")
                        sz = nis_get_size()
                        ibuff = nis_get_integer()

                    else
                    	call print_error_message(-1,"Missing record at NIS (number of states)")
                    end if
                    write(cbuff,'(a,i5)') "number of states -  ", ibuff
					if (ibuff /= current_s_block%N1) then
						call print_error_message(-1,"Number of states in MEF not equal to the number of states in SSF")
					end if

                    !
                    ! Real part of the initial condition
                    !
					if (nis_has_next()) then
                    	call nis_next(err)
                        call print_error_message(err,"reading initial condition from NIS")
                        sz = nis_get_size()

                        allocate(rbuff_v(sz(1)))
                        rbuff_v = nis_get_real(sz(1))

                        do i = 1, sz(1)
                        	ini_dm(i,i) = rbuff_v(i)
						end do

						deallocate(rbuff_v)

                    else
                    	call print_error_message(-1,"Missing record at NIS (initial condition)")
                    end if
                    write(cbuff,'(a,3i7)') "initial condition (Re) read successfully "
                    call print_log_message(trim(cbuff),7)

                    !
                    ! Imaginary part of the initial condition
                    !
					if (nis_has_next()) then
                    	call nis_next(err)
                        call print_error_message(err,"reading initial condition from NIS")
                        sz = nis_get_size()

                        allocate(rbuff_v(sz(1)))
                        rbuff_v = nis_get_real(sz(1))

                        do i = 1, sz(1)
                        	ini_dm(i,i) = ini_dm(i,i) + (0.0d0,1.0d0)*rbuff_v(i)
						end do

						deallocate(rbuff_v)

                    else
                    	call print_error_message(-1,"Missing record at NIS (initial condition)")
                    end if
                    write(cbuff,'(a,3i7)') "initial condition (Re) read successfully "
                    call print_log_message(trim(cbuff),7)


            	do

                	if (nis_has_next()) then
                    	call nis_next(err)
                        call print_error_message(err,"reading NIS at basis type")
                        sz = nis_get_size()
                        write(cbuff,'(i3)') sz(1)
                        write(caux,'('//trim(cbuff)//'a)') nis_get_string(nis_get_size())
                    else
                    	call print_error_message(-1,"Missing record at NIS (basis type")
                    end if
                    write(cbuff,'(a,a)') "basis type -  ", trim(caux)
                    call print_log_message(trim(cbuff),7)

            		if (trim(caux) == "DM_INITIAL_CONDITION_END") exit

				end do

            else if (trim(caux) == "END_SECTIONS") then

                exit

            else

                    call print_error_message(-1,"wrong section separator")

            end if

        end do

        call resources_rewind_blocks()
        call resources_rewind_gofts()
        call resources_set_all_pointers(RESOURCES_SITE_REP)

    end subroutine read_input


end module input_wizard
!
! End of file
!
