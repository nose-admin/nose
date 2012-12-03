!
!  Module handling calculation and loading of relaxation rates
!
!  It provides access to Redfield and Foerster rate theories
!
!
module rates

    use goft
    use sci_foerster
    use sci_redfield

    implicit none

	integer :: load_rates   ! 0 = no, -1 = save, 1 = read
	character(len=1) :: blocknr

contains

    !
    ! initialization of the rates
    !
    subroutine init_rates(err)
        integer, intent(out) :: err
        integer :: i,j,k,jj,kk,a,b, Na, Nb, Nn, fi, Nc
        real(dp) :: sm
        logical :: old

        real(dp), dimension(:,:), pointer :: JJinter
        real(dp), dimension(:,:), allocatable :: Om
        real(dp), dimension(:,:), allocatable :: cw,cw_s

        err = 0

        if ((parallel_id == 0).and.(i_realization == 1)) then

        open(unit=11,file=trim(file_join("log","rates_relax.dat")))
        write(11,'(a)')    "#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        write(11,'(a)')    "# "
        write(11,'(a)')    "#               Energy transfer rates "
        write(11,'(a)')    "# "
        write(11,'(a)')    "#           Calculated by NOSE ver. 0.5.9"
        write(11,'(a)')    "# "
        write(11,'(a)')    "#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        write(11,'(a)')    "# "
        write(11,'(a)')    "# Number of blocks "
        write(11,'(a)')    "# -----------------------------------------------------"
        write(11,'(i5)')    nr_blocks
        write(11,'(a)')    "# "

        end if

		load_rates = 0

        do a = 1, nr_blocks
            do b = 1, nr_blocks

				!
				! Inside blocks we use Redfield theory
				!
                if (a == b) then

                    current_e_block => iblocks(a,a)%eblock
                    call resources_set_all_pointers(RESOURCES_EXCITON_REP)
                    Nn = current_e_block%N1

                    write(blocknr,'(i1)') a

					!
					! frequencies of excitonic transitions
					!
                    allocate(Om(Nn,Nn))
					do i = 1,Nn
						do j = 1,Nn
							Om(i,j) = iblocks(a,a)%eblock%en(i) - iblocks(a,a)%eblock%en(j) &
									-(iblocks(a,a)%eblock%ll(i) - iblocks(a,a)%eblock%ll(j))
						end do
					end do

					! FTs of correlation functions at excitonic frequencies
					allocate(cw(Nn,Nn*Nn))

					! check how many gofts there is
					call resources_rewind_gofts()
					jj = 0
					do
						jj = jj + 1
           		 		if (.not.resources_have_next_goft()) exit
            			call resources_next_goft()
        			end do
					Nc = jj  ! number of defined g(t)s
					allocate(cw_s(Nc,Nn*Nn))
					cw_s = 0.0_dp

					call resources_rewind_gofts()
					jj = 0
      				do

						jj = jj + 1
						! loop over modes
		            	do k = 1, size(current_s_goft%types)

						!	if (current_s_goft%types(k) == "BROWNIAN") then

                        !		do i = 1,Nn
                        !		do j = 1,Nn
                        !			print *, 'O ', comega_overdamped_brownian(Om(i,j), &
                        !			   current_s_goft%params(1,k),current_s_goft%params(2,k),temp*kB_intK)
                        !			print *, 'U ', comega_general_brownian(Om(i,j), &
                        !			   current_s_goft%params(1,k),current_s_goft%params(2,k),0.0d0,temp*kB_intK)
						!		end do
						!		end do

						!	end if

		            	!end do

						!stop

							if (current_s_goft%types(k) == "BROWNIAN") then
                        		!print *, current_s_goft%types(k), current_s_goft%params(1:2,k)

                        		! calculate FT of the correlation function at all excitonic transitions
                        		kk = 0
                        		do i = 1,Nn
                        		do j = 1,Nn
                        			kk = kk + 1
                        			cw_s(jj,kk) = cw_s(jj,kk) + &
                        			comega_overdamped_brownian(Om(i,j), &
                        			   current_s_goft%params(1,k),current_s_goft%params(2,k),temp*kB_intK)
								end do
								end do

							else if (current_s_goft%types(k) == "BROWNIAN_GENERAL") then

                        		kk = 0
                        		do i = 1,Nn
                        		do j = 1,Nn
                        			kk = kk + 1
                        			cw_s(jj,kk) = cw_s(jj,kk) + &
                        			comega_general_brownian(Om(i,j), &
                        			   current_s_goft%params(1,k),current_s_goft%params(2,k),  &
                        			   current_s_goft%params(3,k),temp*kB_intK)

								end do
								end do


							end if

                        end do

           		 		if (.not.resources_have_next_goft()) exit

            			call resources_next_goft()

        			end do

					! give the correlation fucntion values to sites
					do k = 1,Nn
					jj = 0
					do i = 1,Nn
						do j = 1,Nn
							jj = jj + 1

							!print *, iblocks(a,a)%sblock%gindex(k)
							cw(k,jj) = cw_s(iblocks(a,a)%sblock%gindex(k),jj);

						end do
					end do
					end do


!		        	write(cbuff,'(a,i5,a)') "Total of ",N_realizations_local," on each processor"
       				call print_log_message("Calling Redfield rates",6)!trim(cbuff),5)

!       				open(unit=74,file="ahoj.dat")
!
!					do k = 1,Nn
!					jj = 0
!					do i = 1,Nn
!					do j = 1,Nn
!					jj = jj + 1
 !      				write(74,*) cw(k,jj)
  !     				end do
   !    				end do
    !   				end do
     !  				close(74)

					call init_redfield_omega(Om,SS,cw,temp)

! ----------------------------------------------------------------------------------------------------------------------

					if (.false.) then

					if ((load_rates == 0).or.(load_rates == -1)) then

						! rates are calcuted
                    	call init_redfield(en,ll,SS,cc,gt(1)*dt,temp)

						! rates are saved
						if (load_rates == -1) then

							open(unit=74,file=trim(file_join("log","redfield_load.dat")))

							do i = 1,Nn
								write(74,*) redfield_deph_01ex(i)
							end do
							do i = 1,Nn
								do j = 1,Nn
									write(74,*) redfield_deph_1ex(i,j)
								end do
							end do

							do i = 1,Nn
								do j = 1,Nn
									write(74,*) redfield_rates_1ex(i,j)
								end do
							end do

							close(74)

						end if

					else if (load_rates == 1) then

						! rates are loaded
						call init_redfield(en)    ! not calculated, only allocated

						call open_aux_file(74,trim(file_join("log","redfield_load.dat")),"REDFIELD"//blocknr,err)

						read(74,*) fi

						if (fi==1) then

							do i = 1,Nn
								read(74,*) redfield_deph_01ex(i)
							end do
							do i = 1,Nn
								do j = 1,Nn
									read(74,*) redfield_deph_1ex(i,j)
								end do
							end do
							do i = 1,Nn
								do j = 1,Nn
									read(74,*) redfield_rates_1ex(i,j)
								end do
							end do

						else if (fi==-1) then

							redfield_deph_01ex = 0.0d0
							redfield_deph_1ex = 0.0d0
							redfield_rates_1ex = 0.0d0
						end if

	                    call close_aux_file(74)

					end if

					end if

! -----------------------------------------------------------------------------------------------------------



                    nullify(current_e_block%rr)
                    nullify(current_e_block%dr)
                    nullify(current_e_block%dr0)
                    allocate(current_e_block%rr(Nn,Nn))
                    allocate(current_e_block%dr(Nn,Nn))
                    allocate(current_e_block%dr0(Nn))

                    current_e_block%rr = redfield_rates_1ex
                    current_e_block%dr = redfield_deph_1ex
                    current_e_block%dr0 = redfield_deph_01ex

                    if ((parallel_id == 0).and.(i_realization == 1)) then

                        write(11,'(a)')    "#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                        write(11,'(a)')    "# "
                        write(11,'(a)')    "# Relaxation inside a block"
                        write(11,'(a)')    "# "
                        write(11,'(a)')    "# Block id   number of states "
                        write(11,'(a)')    "# -----------------------------------------------------"
                        write(11,'(i7,i15)')    current_e_block%id, current_e_block%N1
                        write(11,'(a)')    "# "
                        write(11,'(a)')    "# Method of calculation: QME (Redfield)"
                        write(11,'(a)')    "#                        2nd order in system-bath cplng"
                        write(11,'(a)')    "#                        Markov approximation          "
                        write(11,'(a)')    "#                        Secular approximation "
                        write(11,'(a)')    "# "
                        write(11,'(a)')    "# Relaxation rates (and corresponding relaxation times): "
                        write(11,'(a)')    "# "
                        write(11,'(a)')    "# state1  state2       rate [1/fs]          time [fs] "
                        write(11,'(a)')    "# -----------------------------------------------------"
                        do j = 1,N1
                            do i = 1, N1

                            if (abs(current_e_block%rr(i,j)) >= 1.0e-10_dp) then
                            write(11,'(i5,a3,i5,a3,e18.6,f18.2)') i,"   ",j,"   ", &
                                current_e_block%rr(i,j), 1/current_e_block%rr(i,j)
                            else
                            write(11,'(i5,a3,i5,a3,e18.6,f18.2)') i,"   ",j,"   ", &
                                current_e_block%rr(i,j), 0.0_dp
                            end if
                            end do
                        end do

                        write(11,'(a)')    "# "
                        write(11,'(a)')    "# Diagonal elements (should be grater then 0): "
                        write(11,'(a)')    "# "
                        write(11,'(a)')    "# state1  state2       rate [1/fs]          time [fs] "
                        write(11,'(a)')    "# -----------------------------------------------------"
                        do i = 1,N1
                            !do j = 1, N1
                            j = i
                            if (abs(current_e_block%rr(i,j)) >= 1.0e-10_dp) then
                            write(11,'(i5,a3,i5,a3,e18.6,f18.2)') i,"   ",j,"   ", &
                                current_e_block%rr(i,j), 1/current_e_block%rr(i,j)
                            else
                            write(11,'(i5,a3,i5,a3,e18.6,f18.2)') i,"   ",j,"   ", &
                                current_e_block%rr(i,j), 0.0_dp
                            end if
                            !end do
                        end do

                        write(11,'(a)') "# "
                        write(11,'(a)') "# "

                        write(11,'(a)') "# "
                        write(11,'(a)') "# Column sums (should be 0.0): "
                        write(11,'(a)') "# "
                        write(11,'(a)') "# column         sum "
                        write(11,'(a)') "# ---------------------------"
                        do j = 1, N1
                            sm = sum(current_e_block%rr(1:N1,j))
                            write(11,'(i5,a3,2e18.6)') j,"   ", sm
                        end do
                        write(11,'(a)') "# "
                        write(11,'(a)') "# "


                    end if


				!
				! Between blocks it is Foerster rates which is used
				!
                else


                    Na = iblocks(a,a)%eblock%N1 ! acceptor block
                    Nb = iblocks(b,b)%eblock%N1 ! donnor block

                    allocate(Om(Na,Nb))
                    do i = 1, Na
                        do j = 1, nb
                           Om(i,j) = iblocks(a,a)%eblock%en(i)-iblocks(b,b)%eblock%en(j)
                        end do
                    end do

                    JJinter => iblocks(a,b)%sblock%J

                    sm = gt(1)*dt
                    call init_foerster(Na,Nb,JJinter, &
                                                transpose(iblocks(a,a)%eblock%gg), &
                                                transpose(iblocks(b,b)%eblock%gg),  &
                                                iblocks(a,a)%eblock%ll,iblocks(b,b)%eblock%ll,  &
                                                Nt(1),sm,Om,Temp)


                    deallocate(Om)

                    if (associated(iblocks(a,b)%sblock%rr)) then
                            iblocks(a,b)%sblock%rr => null() !deallocate(current_e_block%rr)
                    end if
                    if (.not.associated(iblocks(a,b)%sblock%rr)) then
                        allocate(iblocks(a,b)%sblock%rr(Na,Nb))
                    end if


                    if ((parallel_id == 0).and.(i_realization == 1)) then

                        write(11,'(a)')    "#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                        write(11,'(a)')    "# "
                        write(11,'(a)')    "# Relaxation between blocks"
                        write(11,'(a)')    "# "
                        write(11,'(a)')    "# (block ids) (nrs. of states) "
                        write(11,'(a)')    "#   id1  id2     nr1  nr2 "
                        write(11,'(a)')    "# -----------------------------------------------------"
                        write(11,'(i6,i5,i8,i5)') iblocks(a,a)%eblock%id, iblocks(b,b)%eblock%id, &
                                                   iblocks(a,a)%eblock%N1, iblocks(b,b)%eblock%N1
                        write(11,'(a)')    "# "
                        write(11,'(a)')    "# Method of calculation: QME (Foerster)"
                        write(11,'(a)')    "#                        2nd order in resonance cplng"
                        write(11,'(a)')    "#                        Markov approximation          "
                        write(11,'(a)')    "#                        Secular approximation "
                        write(11,'(a)')    "# "
                        write(11,'(a)')    "# Relaxation rates (and corresponding relaxation times): "
                        write(11,'(a)')    "# "
                        write(11,'(a)')    "# state1  state2       rate [1/fs]          time [fs] "
                        write(11,'(a)')    "# -----------------------------------------------------"
                        do i = 1,Na
                            do j = 1, Nb

                            if (abs(foerster_rate(i,j)) >= 1.0e-10_dp) then
                            write(11,'(i5,a3,i5,a3,e18.6,f18.2)') i,"   ",j,"   " , &
                                foerster_rate(i,j), 1/foerster_rate(i,j)
                                !iblocks(1,1)%eblock%rr(i,j), 1/iblocks(1,1)%eblock%rr(i,j)
                            else
                            write(11,'(i5,a3,i5,a3,e18.6,f18.2)') i,"   ",j,"   " , &
                                foerster_rate(i,j), 0.0
                            end if

                            iblocks(a,b)%sblock%rr(i,j) = foerster_rate(i,j)

                            end do
                        end do


                        write(11,'(a)') "# "
                        write(11,'(a)') "# Column sums (should be 0.0): "
                        write(11,'(a)') "# "
                        write(11,'(a)') "# column         sum "
                        write(11,'(a)') "# ---------------------------"
                        sm = 0.0_dp
                        do j = 1, Nb
                            sm = sum(iblocks(a,b)%sblock%rr(1:Na,j))
                            write(11,'(i5,a3,2e18.6)') j,"   " , sm
                        end do
                        write(11,'(a)') "# "
                        write(11,'(a)') "# "

                        call clean_foerster()

                    end if

                end if

            end do
        end do

        if ((parallel_id == 0).and.(i_realization == 1)) then
            write(11,'(a)')    "# END OF OUTPUT "
        end if

        close(unit=11)

    end subroutine init_rates


    subroutine clean_rates



    end subroutine clean_rates

end module rates
