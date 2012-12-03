!
! Module TDPT-3
!

module module_tdpt3

#include <util_allocation.h>
    use prepare
    use twod
    use pump_probe
    use numer_fft

    implicit none

    public  :: do_tdpt3_work

    ! time dependent signals
    complex, dimension(:), allocatable  :: polar_1
    complex, dimension(:), allocatable  :: polar_fl
    complex, dimension(:), allocatable  :: tdep_cd

    ! spectra
    real, dimension(:), allocatable     :: spect_abs
    real, dimension(:), allocatable     :: spect_fluor
    real, dimension(:), allocatable     :: spect_cd

    ! spectra
    complex, dimension(:,:,:), allocatable :: spect_2D_ftpe
    complex, dimension(:,:,:), allocatable :: spect_2D_gs
    complex, dimension(:,:,:), allocatable :: spect_2D_esa


    ! statistics
    real, dimension(:), allocatable     :: hist_relax
    integer, parameter :: N_hist = 300

	! populations
	real, dimension(:,:), allocatable        :: ppop


    ! Fourier transform
    integer             :: NFFT, padfac
    integer, parameter :: NFFT_basic = 1024


    ! privates
    private  :: dom
    real(dp) :: dom

    real(dp), private ::    hom_l
    real(dp), private ::    hom_u
    real(dp), private ::    dhom

    character(len=256), private :: char_buff

    character(len=10), private :: sbuff, sbuff2

contains


    !
    ! Do all the simulations within the module
    !
    subroutine do_tdpt3_work(err)
        integer, intent(out) :: err

        integer :: i
        real(dp) :: t

        logical :: have_polar_1, have_polar_fl, have_2d_ftpe, have_pump_probe

        have_polar_1 = .false.
        have_polar_fl = .false.
        have_2d_ftpe = .false.
        have_pump_probe = .false.

        err = 0

        call prepare_bath_interactions(err)
        call update_resources_tdpt3()

        call resources_rewind_blocks()
        call resources_set_all_pointers(RESOURCES_EXCITON_REP)
        call rewind_evolution_operators()

        ! calculate first order polarization

        Ueg => current_Ueg%Ueg

        !*************************************************************
        !  Calculation of the first order polarization
        !*************************************************************

        if (resources_output_contains("polar_1")) then

            ! output of the first order polarization

            call create_polar_1()
            have_polar_1 = .true.


        end if

        !*************************************************************
        !  Calculation of the first order polarization
        !*************************************************************

        if (resources_output_contains("polar_fl")) then

            ! output of the first order polarization

            call create_polar_fl()
            have_polar_fl = .true.


        end if


        !*************************************************************
        ! Calculation of linear absorption spectrum
        !*************************************************************
        if (resources_output_contains("spect_abs")) then

            ! output linear absorption spectrum

            if (.not.have_polar_1) then

                call create_polar_1()
	            have_polar_1 = .true.

            end if


        end if


        !*************************************************************
        ! Calculation of the circular dichroism spectrum
        !*************************************************************

        if (resources_output_contains("spect_CD")) then

            ! create circular dichroism spectrum
            call create_tdep_CD()

        end if



        if (resources_output_contains("spect_fluor")) then

            if (.not.have_polar_fl) then

                call create_polar_fl()

            end if

        end if

        if (resources_output_contains("hist_relax")) then

            call create_hist_relax()

        end if

        !*************************************************************
        !  Calculation of 2D spectrum
        !*************************************************************
        if (resources_output_contains("2d_ftpe")) then

            ! output of populations and coherences

            call create_2D_spec()
            have_2d_ftpe = .true.


        end if


        !*************************************************************
        !  Calculation of pump-probe spectrum
        !*************************************************************
        if (resources_output_contains("pump_probe")) then

            ! output of populations and coherences

            call create_pump_probe()
            have_pump_probe = .true.


            !stop


        end if

        !*************************************************************
        ! Calculation of population dynamics
        !*************************************************************
        if (resources_output_contains("dens_exc_block")) then

            ! output linear absorption spectrum
            if ((parallel_id == 0).and.(i_realization == 1)) then
        		call create_populations()
			end if

        end if



       ! call resources_clean_excitons()
       ! call clean_resources_tdpt3()



    end subroutine do_tdpt3_work



    !
    ! Cleans all allocated memory
    !
    subroutine clean_tdpt3

    end subroutine clean_tdpt3

    !
    ! Collecting data after simulation finished
    !
    subroutine collect_tdpt3_data(err)
        integer, intent(out) :: err
        integer :: i,j
        real(dp) :: oma, rmx

        err = 0

        !*************************************************************
        !  Output of the first order polarization
        !*************************************************************

        if (resources_output_contains("polar_1")) then

            ! output of the first order polarization

            call collect_polar_1()


            if (parallel_id == 0) then

            	write(char_buff,*) i_realization*parallel_nr
				sbuff2 = trim(adjustl(char_buff))

                open(unit=11,file=trim(file_join(out_dir,"polar_1.dat_"//trim(sbuff2))))

                 do i = 1, Nt(1)

                    write(11,*) i*gt(1)*dt, real(polar_1(i)), aimag(polar_1(i))

                 end do

                close(unit=11)
            end if

			if (.not.resources_output_contains("spect_abs")) then
            	polar_1 = 0.0_dp
			end if

        end if


        if (resources_output_contains("polar_fl")) then

            ! output of the first order polarization

            call collect_polar_fl()

            if (parallel_id == 0) then
                open(unit=11,file=trim(file_join(out_dir,"polar_fl.dat")))

                 do i = 1, Nt(1)

                    write(11,*) i*gt(1)*dt, real(polar_fl(i)), aimag(polar_fl(i))

                 end do

                close(unit=11)
            end if

        end if


        !*************************************************************
        ! Output of linear absorption spectrum
        !*************************************************************
        if (resources_output_contains("spect_abs")) then

            ! output linear absorption spectrum

            call create_spect_abs()

            if (parallel_id == 0) then

            	write(char_buff,*) i_realization*parallel_nr
				sbuff2 = trim(adjustl(char_buff))

                open(unit=11,file=trim(file_join(out_dir,"spect_abs.dat_"//trim(sbuff2))))

				if (cur_outnr > 0 ) then
                	!print *, "Outputs: ", trim(outppars(cur_outnr,1))
                	if (trim(outppars(cur_outnr,1))=="normalized") then
                		rmx = maxval(spect_abs)
                		spect_abs = spect_abs/rmx
						write(char_buff,'(a)') "   Normalizing output "
    					call print_log_message(trim(char_buff),5)
                	end if
				end if

                do i = 1, NFFT
                    oma = ((-NFFT/2 + i)*dom + rwa)*Energy_internal_to_cm

                    write(11,*) oma, spect_abs(i)

                end do


                close(unit=11)
            end if

			spect_abs = 0.0_dp
			polar_1 = 0.0_dp

        end if

        !*************************************************************
        ! Output of linear absorption spectrum
        !*************************************************************
        if (resources_output_contains("spect_fluor")) then

            ! output linear absorption spectrum

            call create_spect_fluor()

            if (parallel_id == 0) then
                open(unit=11,file=trim(file_join(out_dir,"spect_fluor.dat")))

				if (cur_outnr > 0 ) then
                	!print *, "Outputs: ", trim(outppars(cur_outnr,1))
                	if (trim(outppars(cur_outnr,1))=="normalized") then
                		rmx = maxval(spect_fluor)
                		spect_fluor = spect_fluor/rmx
						write(char_buff,'(a)') "   Normalizing output "
    					call print_log_message(trim(char_buff),5)
                	end if
				end if

                do i = 1, NFFT
                    oma = ((-NFFT/2 + i)*dom + rwa)*Energy_internal_to_cm

                    write(11,*) oma, spect_fluor(i)

                end do


                close(unit=11)
            end if

        end if

        !*************************************************************
        ! Output of the circular dichroism spectrum
        !*************************************************************


        if (resources_output_contains("spect_CD")) then

            ! output circular dichroism spectrum

            call create_spect_CD()

            if (parallel_id == 0) then
                open(unit=11,file=trim(file_join(out_dir,"spect_CD.dat")))

				if (cur_outnr > 0 ) then
                	!print *, "Outputs: ", trim(outppars(cur_outnr,1))
                	if (trim(outppars(cur_outnr,1))=="normalized") then
                		rmx = maxval(spect_CD)
                		if (rmx == 0.0_dp) then
                			spect_CD = 0.0_dp
                	    else
                			spect_CD = spect_CD/rmx
							write(char_buff,'(a)') "   Normalizing output "
    						call print_log_message(trim(char_buff),5)
                		end if
                	end if
				end if

                do i = 1, NFFT
                    oma = ((-NFFT/2 + i)*dom + rwa)*Energy_internal_to_cm

                    write(11,*) oma, spect_cd(i)

                end do


                close(unit=11)
            end if

        end if

        !****************************************************************
        ! Output of relaxation rate spectrum
        !****************************************************************
        if (resources_output_contains("hist_relax")) then

            ! output circular dichroism spectrum

            call collect_hist_relax()

            if (parallel_id == 0) then
                open(unit=11,file=trim(file_join(out_dir,"hist_relax.dat")))

                do i = 1, N_hist
                    !oma = (-(-NFFT/2 + i)*dom + rwa)*Energy_internal_to_cm

                    write(11,*) (hom_l + i*dhom)*Energy_internal_to_cm, hist_relax(i)

                end do


                close(unit=11)
            end if

        end if


        !*************************************************************
        !  Outputting 2D spectra
        !*************************************************************
        if (resources_output_contains("2d_ftpe")) then

            call collect_2D_spec()

            rmx = maxval(real(spect_2d_ftpe(:,:,1)))
            spect_2d_ftpe = spect_2d_ftpe/rmx
            spect_2d_esa = spect_2d_esa/rmx
			spect_2d_gs = spect_2d_gs/rmx

            dom = 2.0_dp*PI_D/(NFFT*gt(1)*dt)

            if (parallel_id == 0) then

            	do i = 1, Nt(2)

				write(char_buff,*) nint(dt*gt(2)*(i-1))
				sbuff = trim(adjustl(char_buff))
				write(char_buff,*) i_realization*parallel_nr
				sbuff2 = trim(adjustl(char_buff))

                call save_2D_tot(real(spect_2d_ftpe(:,:,i)),trim(file_join(out_dir,"twod_re_tot_T=" &
                        //trim(sbuff)//"fs.dat_"//trim(sbuff2))))
                call save_2D_tot(aimag(spect_2d_ftpe(:,:,i)),trim(file_join(out_dir,"twod_im_tot_T=" &
                		//trim(sbuff)//"fs.dat_"//trim(sbuff2))))
                call save_2D_tot(real(spect_2d_esa(:,:,i)),trim(file_join(out_dir,"twod_re_esa_T=" &
                		//trim(sbuff)//"fs.dat_"//trim(sbuff2))))
                call save_2D_tot(aimag(spect_2d_esa(:,:,i)),trim(file_join(out_dir,"twod_im_esa_T=" &
                		//trim(sbuff)//"fs.dat_"//trim(sbuff2))))
                call save_2D_tot(real(spect_2d_gs(:,:,i)),trim(file_join(out_dir,"twod_re_gs_T="   &
                		//trim(sbuff)//"fs.dat_"//trim(sbuff2))))
                call save_2D_tot(aimag(spect_2d_gs(:,:,i)),trim(file_join(out_dir,"twod_im_gs_T="   &
                		//trim(sbuff)//"fs.dat_"//trim(sbuff2))))
                call save_2D_limits(dom,NFFT,trim(file_join(out_dir,"twod_limits.dat")))

				end do

            end if

            spect_2d_ftpe = 0.0_dp
            spect_2d_esa = 0.0_dp
			spect_2d_gs = 0.0_dp


        end if



        !*************************************************************
        !  Outputting pump-probe spectra
        !*************************************************************
        if (resources_output_contains("pump_probe")) then

            call collect_pump_probe()

            dom = 2.0_dp*PI_D/(NFFT*gt(1)*dt)

            if (parallel_id == 0) then

              !  call save_2D_tot(real(spect_2d_ftpe),trim(file_join(out_dir,"twod_re_tot_T=0fs.dat")))
              !  call save_2D_tot(aimag(spect_2d_ftpe),trim(file_join(out_dir,"twod_im_tot_T=0fs.dat")))
              !  call save_2D_limits(dom,NFFT,trim(file_join(out_dir,"twod_limits.dat")))

            end if

        end if

        if (resources_output_contains("dens_exc_block")) then

            ! output of the first order polarization

            if ((parallel_id == 0)) then
                open(unit=11,file=trim(file_join(out_dir,"dens_exc_block.dat")))

                do i = 1, Nt(2)

                    write(11,*) (i-1)*gt(2)*dt, ppop(:,i)

                end do

				!do i = 1, N1
				!do j = 1, N1
				!
				!	write(11,*) i,j, evops(1,1)%Uee(i,i,j,j,1)
				!
				!end do
				!end do

                close(unit=11)

                DEALLOCATE(ppop)

            end if



        end if



    end subroutine collect_tdpt3_data


    !
    !
    !


    !
    ! Calculation of the first order polarization
    !
    subroutine create_polar_1()
        integer(i4b) :: i,k,b, bb
		character(len=64) :: cbuff

        if (.not.allocated(polar_1)) then
            allocate(polar_1(1:Nt(1)))
            polar_1 = 0.0_dp
        end if

        call resources_rewind_blocks()
        call resources_set_all_pointers(RESOURCES_EXCITON_REP)
        call rewind_evolution_operators()


        write(cbuff,'(a)') "Creating linear polarizations"
        call print_log_message(trim(cbuff),6)

        ! loop over blocks
        bb = 0
        do

        	bb = bb + 1
	        Ueg => current_Ueg%Ueg

			do k = 1, N1
            	polar_1(1) = polar_1(1) + (abs(current_e_block%dd(k))**2)
            end do
            ! loop over time
            do i = 1, Nt(1)-1
                ! loop over excited states
                do k = 1, N1
                    polar_1(i+1) = polar_1(i+1) + (abs(current_e_block%dd(k))**2)*Ueg(k,i+1)
                end do
            end do

			if (maxval(abs(polar_1)) > abs(polar_1(1))) then


            	write(char_buff,'(a,i5,a)') "Broken realization nr. ", &
              	i_realization, ". Thrown away ..."
            	call print_warning_message(trim(char_buff),5)
            	polar_1 = 0.0_dp

			end if

            if (.not.resources_have_next_block()) exit

            call resources_next_block()
            call resources_set_all_pointers(RESOURCES_EXCITON_REP)
            call next_evolution_operators()

        end do

    end subroutine create_polar_1

    !
    ! Calculation of the polarization for fluorescence
    !
    subroutine create_polar_fl()
        integer(i4b) :: i,k,b

        if (.not.allocated(polar_fl)) then
            allocate(polar_fl(1:Nt(1)))
            polar_fl = 0.0_dp
        end if

        call resources_rewind_blocks()
        call resources_set_all_pointers(RESOURCES_EXCITON_REP)
        call rewind_evolution_operators()
        Ueg => current_Ueg%feg


        ! loop over blocks
        do

            ! loop over time
            do i = 1, Nt(1)
                ! loop over excited states
                do k = 1, N1
                    polar_fl(i) = polar_fl(i) + current_e_block%fl_fac(k)*(abs(current_e_block%dd(k))**2)*Ueg(k,i)
                end do
            end do

            if (.not.resources_have_next_block()) exit

            call resources_next_block()
            call resources_set_all_pointers(RESOURCES_EXCITON_REP)
            call next_evolution_operators()
            Ueg => current_Ueg%feg

        end do

    end subroutine create_polar_fl


    !
    !
    subroutine collect_polar_1()
        integer :: N
        character(len=256) :: buff

        N = N_since_saved !N_realizations_local
        call parallel_average_complex_array(polar_1,N)
        polar_1_collected = .true.

        write(buff,'(a,i5,a)') "Polar_1  : Total of ",N," realization(s) received"
        call print_log_message(trim(buff),5)

    end subroutine collect_polar_1

    !
    !
    subroutine collect_polar_fl()
        integer :: N
        character(len=256) :: buff

        N = N_realizations_local
        call parallel_average_complex_array(polar_fl,N)
        polar_fl_collected = .true.

        write(buff,'(a,i5,a)') "Polar_fl : Total of ",N," realization(s) received"
        call print_log_message(trim(buff),5)

    end subroutine collect_polar_fl


    !
    ! Calculate linear absorption spectrum from polarization
    !
    subroutine create_spect_abs()

        complex(dpc), dimension(:),   allocatable    :: sig
        complex(dpc), dimension(:,:), allocatable    :: dat
        real(dp), dimension(:), allocatable :: rs

        integer :: on, ch, i
        real(dp) :: oma, rr

        padfac = 4
        NFFT = (2**padfac)*NFFT_basic

        allocate(sig(NFFT))
        allocate(dat(1,NFFT))
        allocate(rs(NFFT))

       	if (.not.allocated(spect_abs)) then
        	allocate(spect_abs(NFFT))
        end if

        rs        = 0.0_dp
        spect_abs = 0.0_dp

        dom = 2.0_dp*PI_D/(NFFT*gt(1)*dt)
        !print *, "***************** ", NFFT, gt(1), dt, dom, rwa

        if (.not.polar_1_collected) then
            call collect_polar_1()
        end if


        sig = 0.0_dp
        sig(1:Nt(1)) = polar_1(1:Nt(1))

        do i = 1, NFFT/2
            sig(NFFT/2+i) = conjg(sig(NFFT/2-i+2))
        end do

        !
        ! FFT a la numerical recipes
        !
        dat(1,:) = sig(:)
        call fft_row(dat,1_i4b)
        sig(1:NFFT/2-1)  = dat(1,NFFT/2+2:NFFT)
        sig(NFFT/2:NFFT) = dat(1,1:NFFT/2+1)


        do i = 1, NFFT
            oma = (-NFFT/2 + i)*dom + rwa

!                if (i == 0) then
!                    rs(i) = abs(oma)*real(sig(NFFT))
!                else
                    rs(i) = abs(oma)*real(sig(i))
!                end if


        end do

        !rr = maxval(rs(:))

        ! no normalization
        !rr = 1.0_dp

        spect_abs = rs !/rr

        ! we normalize at the output

    end subroutine create_spect_abs

    !
    !
    subroutine collect_spect_abs()


    end subroutine collect_spect_abs


    !
    ! Calculate fluorescence spectrum from polarization
    !
    subroutine create_spect_fluor()

        complex(dpc), dimension(:),   allocatable    :: sig
        complex(dpc), dimension(:,:), allocatable    :: dat
        real(dp), dimension(:), allocatable :: rs

        integer :: on, ch, i
        real(dp) :: oma, rr

        padfac = 4
        NFFT = (2**padfac)*NFFT_basic

        allocate(sig(NFFT))
        allocate(dat(1,NFFT))
        allocate(rs(NFFT))
        allocate(spect_fluor(NFFT))
        rs        = 0.0_dp
        spect_fluor = 0.0_dp

        dom = 2.0_dp*PI_D/(NFFT*gt(1)*dt)
        !print *, NFFT, dt, dom

        if (.not.polar_fl_collected) then
            !print *, "CC"
            call collect_polar_fl()
        end if


        sig = 0.0_dp
        sig(1:Nt(1)) = polar_fl(1:Nt(1))

        do i = 1, NFFT/2
            sig(NFFT/2+i) = conjg(sig(NFFT/2-i+2))
        end do

        !
        ! FFT a la numerical recipes
        !
        dat(1,:) = sig(:)
        call fft_row(dat,1_i4b)
        sig(1:NFFT/2-1)  = dat(1,NFFT/2+2:NFFT)
        sig(NFFT/2:NFFT) = dat(1,1:NFFT/2+1)


        do i = 1, NFFT
            oma = (-NFFT/2 + i)*dom + rwa

!                if (i == 0) then
!                    rs(i) = abs(oma)*real(sig(NFFT))
!                else
                    rs(i) = abs(oma)*real(sig(i))
!                end if


        end do

        !rr = maxval(rs(:))

        ! no normalization
        !rr = 1.0_dp

        spect_fluor = rs !/rr

		! we normalize at the output

    end subroutine create_spect_fluor



    subroutine create_tdep_CD()
        integer(i4b) :: i,k,b

        if (.not.allocated(tdep_CD)) then
            allocate(tdep_CD(1:Nt(1)))
            tdep_CD = 0.0_dp
        end if

        call resources_rewind_blocks()
        call resources_set_all_pointers(RESOURCES_EXCITON_REP)
        call rewind_evolution_operators()
        Ueg => current_Ueg%Ueg

        ! loop over blocks
        do

            ! loop over time
            do i = 1, Nt(1)
                ! loop over excited states
                do k = 1, N1
                    tdep_CD(i) = tdep_CD(i) + current_e_block%rm(k)*Ueg(k,i)
                end do
            end do

            if (.not.resources_have_next_block()) exit

            call resources_next_block()
            call resources_set_all_pointers(RESOURCES_EXCITON_REP)
            call next_evolution_operators()
            Ueg => current_Ueg%Ueg

        end do

    end subroutine create_tdep_CD

    !
    !
    subroutine collect_tdep_CD()
        integer :: N
        character(len=256) :: buff

        N = N_realizations_local
        call parallel_average_complex_array(tdep_cd,N)
        tdep_cd_collected = .true.

        write(buff,'(a,i5,a)') "Tdep CD  : Total of ",N," realization(s) received"
        call print_log_message(trim(buff),5)

    end subroutine collect_tdep_CD

    !
    ! Calculate linear absorption spectrum from polarization
    !
    subroutine create_spect_CD()

        complex(dpc), dimension(:),   allocatable    :: sig
        complex(dpc), dimension(:,:), allocatable    :: dat
        real(dp), dimension(:), allocatable :: rs

        integer :: on, ch, i
        real(dp) :: oma, rr

        padfac = 4
        NFFT = (2**padfac)*NFFT_basic

        allocate(sig(NFFT))
        allocate(dat(1,NFFT))
        allocate(rs(NFFT))
        allocate(spect_cd(NFFT))
        rs        = 0.0_dp
        spect_cd  = 0.0_dp

        dom = 2.0_dp*PI_D/(NFFT*gt(1)*dt)
        !print *, NFFT, dt, dom

        if (.not.tdep_cd_collected) then
            call collect_tdep_cd()
        end if


        sig = 0.0_dp
        sig(1:Nt(1)) = tdep_cd(1:Nt(1))

        do i = 1, NFFT/2
            sig(NFFT/2+i) = conjg(sig(NFFT/2-i+2))
        end do

        !
        ! FFT a la numerical recipes
        !
        dat(1,:) = sig(:)
        call fft_row(dat,1_i4b)
        sig(1:NFFT/2-1)  = dat(1,NFFT/2+2:NFFT)
        sig(NFFT/2:NFFT) = dat(1,1:NFFT/2+1)


        do i = 1, NFFT
            oma = (-NFFT/2 + i)*dom + rwa

!                if (i == 0) then
!                    rs(i) = abs(oma)*real(sig(NFFT))
!                else
                    rs(i) = abs(oma)*real(sig(i))
!                end if


        end do

        !rr = maxval(abs(rs(:)))

        ! no normalization
        rr = 1.0_dp

        if (rr /= 0.0_dp) then

            spect_cd = rs/rr

        else

            spect_cd = 0.0

        end if

    end subroutine create_spect_CD

    !
    !
    !
    subroutine collect_spec_CD()

    end subroutine collect_spec_CD

    !
    ! create a local histogram of rates
    !
    subroutine create_hist_relax()

        integer :: k, l, Nhp
        real(dp) :: rr, om

        if (.not.allocated(hist_relax)) then
        allocate(hist_relax(N_hist))
        hist_relax = 0.0_dp
        end if
        hom_l = -2000.0_dp/Energy_internal_to_cm
        hom_u = 2000.0_dp/Energy_internal_to_cm

        dhom = ((hom_u - hom_l)/real(N_hist))

        call resources_rewind_blocks()
        call resources_set_all_pointers(RESOURCES_EXCITON_REP)

        ! loop over blocks
        do

            ! loop over time
            !do i = 1, Nt(1)
                ! loop over excited states
                do k = 1, N1
                    do l = 1, N1
                        if (k /= l) then
                        rr = current_e_block%rr(k,l)
                        om = current_e_block%en(k) - current_e_block%en(l)
                        Nhp = int((om - hom_l)/dhom) + 1
                        if ((Nhp < N_hist).and.(Nhp > 0)) then
                          hist_relax(Nhp) = hist_relax(Nhp) + abs(rr)
                        end if
                        end if
                    end do
                end do
            !end do

            if (.not.resources_have_next_block()) exit

            call resources_next_block()
            call resources_set_all_pointers(RESOURCES_EXCITON_REP)

        end do

    end subroutine create_hist_relax



    subroutine collect_hist_relax()
        integer :: N
        character(len=256) :: buff

        N = N_realizations_local
        call parallel_average_real_array(hist_relax,N)
        hist_relax_collected = .true.

        write(buff,'(a,i5,a)') "Histogram of rates  : Total of ",N," realization(s) received"
        call print_log_message(trim(buff),5)

    end subroutine collect_hist_relax

    !
    !  Calculates 2D spectrum
    !
    subroutine create_2D_spec()

        !if (parallel_id == 0) then
        integer :: Nt1, Nt2
        integer :: kk


        ! create coherences and population evolution

        call create_coherences_eg()
        call create_coherences_fe()
!        call create_populations()
        call create_coherences_ee()


		! loop over population time
		do kk = 1, Nt(2)

		use_Uee_index = kk

        call init_twod()

        Nt1 = size(SR0,1)
        Nt2 = size(SR0,2)

        if (.not.allocated(spect_2d_ftpe)) then
            allocate(spect_2d_ftpe(Nt1,Nt2,Nt(2))) !(size(SR0,1),size(SR0,2)))
            allocate(spect_2d_gs(Nt1,Nt2,Nt(2)))
            allocate(spect_2d_esa(Nt1,Nt2,Nt(2)))
            spect_2d_ftpe = 0.0_dp
            spect_2d_gs = 0.0_dp
            spect_2d_esa = 0.0_dp
        end if

        spect_2d_gs(:,:,kk) = spect_2d_gs(:,:,kk) + SR0 + SNR0
		spect_2d_esa(:,:,kk) = spect_2d_esa(:,:,kk) + ESAR + ESAN
		spect_2d_ftpe(:,:,kk) = spect_2d_gs(:,:,kk) + spect_2d_esa(:,:,kk)

        call clean_twod()

		end do

        !end of the loop over population time

    end subroutine create_2D_spec


    !
    !  Calculates pump-probe spectrum
    !
    subroutine create_pump_probe()

        !if (parallel_id == 0) then
        integer :: Nt1, Nt2

        call create_coherences_eg()
        call create_coherences_fe()
        call create_populations()
        call create_coherences_ee()

        call init_pump_probe()

    end subroutine create_pump_probe


    !
    ! Collects 2D spectra
    !
    subroutine collect_2D_spec()
        integer :: N
        character(len=256) :: buff

        N = N_since_saved ! N_realizations_local
        call parallel_average_complex_array3(spect_2d_esa,N)
        N = N_since_saved !N_realizations_local
        call parallel_average_complex_array3(spect_2d_gs,N)
        N = N_since_saved ! N_realizations_local
        call parallel_average_complex_array3(spect_2d_ftpe,N)


!        spect_2d_ftpe = spect_2d_ftpe/N
        write(buff,'(a,i5,a)') "2d_ftpe  : Total of ",N," realization(s) received"
        call print_log_message(trim(buff),5)

    end subroutine collect_2D_spec


    !
    ! Collects 2D spectra
    !
    subroutine collect_pump_probe()
        integer :: N
        character(len=256) :: buff

        N = N_realizations_local

        !spect_2d_ftpe = spect_2d_ftpe/N
        write(buff,'(a,i5,a)') "pump_probe  : Total of ",N," realization(s) received"
        call print_log_message(trim(buff),5)

    end subroutine collect_pump_probe


    !
    ! Evolution of the optical coherences
    !
    subroutine create_coherences_eg()
    	integer :: k, i

        call resources_rewind_blocks()
        call resources_set_all_pointers(RESOURCES_EXCITON_REP)
        call rewind_evolution_operators()

		!allocate(Uegeg(N1,1,N1,1,Nt(1)))

        ! loop over time
        do i = 1, Nt(1)
       		! loop over excited states
            do k = 1, N1
               evops(1,1)%Ueg(k,1,k,1,i) = current_Ueg%Ueg(k,i)
            end do
        end do


    end subroutine create_coherences_eg

    !
    ! Evolution of the optical coherences
    !
    subroutine create_coherences_fe()


    end subroutine

    !
    ! Evolution of populations
    !
    subroutine create_populations()
		integer :: k,l,t
		real(dp) :: ss
		real(dp), dimension(:), allocatable :: pini

		integer :: pop_ini


		ALLOCATE(ppop,(N1,Nt(2)))
		ALLOCATE(pini,(N1))

		pop_ini = ini_basis_type

!		pop_ini = INI_POPULATION_EXCITON !INI_ULTRA_FAST_PULSE

		if (pop_ini == INI_ULTRA_FAST_PULSE) then

			! normalized initial condition after ultrafast excitation
			ss = 0.0_dp
			do k = 1,N1
				pini(k) = iblocks(1,1)%eblock%dd(k)**2
				ss = ss + pini(k)
			end do
			pini = pini/ss

		else if (pop_ini == INI_POPULATION_EXCITON) then

			ss = 0.0_dp
			do k = 1,N1
				pini(k) = ini_dm(k,k)
				ss = ss + pini(k)
			end do
			pini = pini/ss

		else if (pop_ini == INI_POPULATION_SITES) then

			! transform ini_dm to exciton basis

			ss = 0.0_dp
			do k = 1,N1
				pini(k) = ini_dm(k,k)
				ss = ss + pini(k)
			end do
			pini = pini/ss


		else if (pop_ini == INI_DENSITY_MATRIX) then

		end if

		close(12)

		ppop = 0.0_dp

		! loop over population time
		do t = 1, Nt(2)

			do k = 1, N1
				do l = 1, N1
					ppop(k,t) = ppop(k,t) + evops(1,1)%Uee(k,k,l,l,t)*pini(l)
				end do
			end do

		end do

		DEALLOCATE(pini)

    end subroutine create_populations


    !
    ! Evolution of intraband coherences
    !
    subroutine create_coherences_ee()

    end subroutine



end module module_tdpt3
