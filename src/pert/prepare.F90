!
! Module that takes care of all preparatory
! calculations  relaxation rates, correlation
!                functions etc.
!
! author: Tomas Mancal
! e-mail: tomas.mancal@mff.cuni.cz
!
! last changed: 2006-12-22
!
module prepare

    use bath_interactions
    use orfact

    implicit none

    character(len=1) yn


contains

    subroutine prepare_for_all()

        integer :: err, i,j
        real(dp), dimension(:), allocatable :: ll

        call init_goft(err)

        call resources_rewind_blocks()
        call resources_rewind_gofts()

        allocate(ll(nr_gofts))
        allocate(igofts(nr_gofts))

		i = 1
        do
            ll(current_s_goft%index) = current_s_goft%lambda
            igofts(i)%goft => current_s_goft

            if (.not.resources_have_next_goft()) exit
            call resources_next_goft()

            i = i + 1

        end do

        do

            call resources_set_all_pointers(RESOURCES_SITE_REP)

            ! set values (reorganization energy etc.)
            allocate(current_s_block%ll(N1))
            allocate(current_s_block%dd(N1))

            !print *, "N1 = ", N1

            do i = 1, N1
                ! set ll
               ! print *, "current_s_block%gindex(i) ", i, current_s_block%gindex(i)
                current_s_block%ll(i) = ll(current_s_block%gindex(i))
                current_s_block%dd(i) = sqrt(current_s_block%dx(i)**2 + &
                                             current_s_block%dy(i)**2 + &
                                             current_s_block%dz(i)**2)
            end do

            if (.not.resources_have_next_block()) exit
            call resources_next_block()

        end do

        call resources_rewind_blocks()

        call init_excitons(err)
        call init_orfact()

        ! set block storage
        allocate(iblocks(nr_blocks,nr_blocks))
        allocate(evops(nr_blocks,nr_blocks))


        ! diagonal blocks
        call resources_rewind_blocks()
        i = 1
        do

            iblocks(i,i)%sblock => current_s_block
            iblocks(i,i)%eblock => current_e_block

            if (.not.resources_have_next_block()) exit
            call resources_next_block()
            i = i + 1

        end do

        ! coupling blocks
        call resources_rewind_inter_blocks()
        do i = 1, nr_blocks
            do j = 1, nr_blocks
                if (i /= j) then
                    iblocks(i,j)%sblock => current_i_block
                    if (resources_have_next_inter_block()) then
                        call resources_next_inter_block()
                    end if
                end if
            end do
        end do

        call print_log_message('Checking coupling between blocks',7)
        do i = 1, nr_blocks
            do j = 1, nr_blocks
                if (i /= j) then

                    ! check couplings (fill zeros)
                    if (.not.associated(iblocks(i,j)%sblock%J)) then
                        if (associated(iblocks(j,i)%sblock%J)) then
                            allocate(iblocks(i,j)%sblock%J(size(iblocks(j,i)%sblock%J,1), &
                                                             size(iblocks(j,i)%sblock%J,1)))
                            iblocks(i,j)%sblock%J = transpose(iblocks(j,j)%sblock%J)
                        else
                            allocate(iblocks(i,j)%sblock%J(iblocks(i,i)%sblock%N1, &
                                                             iblocks(j,j)%sblock%N1))
                            iblocks(i,j)%sblock%J = 0.0_dp
                        end if
                    end if

                    ! check if it is symmetric
                    if (maxval(abs(iblocks(i,j)%sblock%J - transpose(iblocks(j,i)%sblock%J)))<=0.0001) then
                        !print *, 'OK'
                    else if (maxval(abs(iblocks(i,j)%sblock%J - 0.0_dp))<=0.0001) then
                        iblocks(i,j)%sblock%J = transpose(iblocks(j,i)%sblock%J)
                        call print_warning_message('Interblock coupling seems to be asymetric',5)
                    else if (maxval(abs(iblocks(j,i)%sblock%J - 0.0_dp))<=0.0001) then
                        iblocks(j,i)%sblock%J = transpose(iblocks(i,j)%sblock%J)
                        call print_warning_message('Interblock coupling seems to be asymetric',5)
                    else
                        call print_error_message(1,'Asymmetry in couplings between blocks')
                        print *, 'Error: '
                        print *, ' '
                        print *, i,j
                        print *, iblocks(i,j)%sblock%J
                        print *, j,i
                        print *, iblocks(j,i)%sblock%J
                    end if

                end if
            end do
        end do
        call print_log_message('...OK',7)

        call init_resources_tdpt3()

    end subroutine prepare_for_all


    !
    ! Call preparation routines in correct order
    !
    ! Here the quantities are prepared for individual realizations
    !
    subroutine prepare_for_one(err)
        integer, intent(out) :: err

        integer :: i,j

        err = 0

        call update_excitons(err)
        call init_orfact()


    end subroutine prepare_for_one


end module prepare
