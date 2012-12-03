#include "util_allocation.h"
!
!
!
!
module module_qme

	use prepare
	use twod

	use resources_qme
	use std_types

	use numer_ode
	use numer_fft
	use sci_misc

	use qme_nakajima_zwanzig
	use qme_weak_excitons

	implicit none

	logical  :: have_populations, have_coherences, have_polar_1, have_pop_coh, &

				 have_2d_ftpe, have_1to2ex_coh, qme_polar_1_collected
	integer  :: NFFT, padfac
	integer, parameter :: NFFT_basic = 1024
	real(dp) :: dom


	! privates
	private :: dom


contains
	!
	! Do all the simulations within the module
	!
	subroutine do_qme_work(err)
		integer, intent(out) :: err

		integer :: i

		have_populations = .false.
		have_coherences = .false.
		have_polar_1 = .false.
		have_pop_coh = .false.
		have_2d_ftpe = .false.
		have_1to2ex_coh = .false.
		qme_polar_1_collected = .false.

		if(.not. qme_module_already_initiated) then
			call init_resources_qme()
			qme_module_already_initiated = .true.
		end if

		call resources_rewind_blocks()  ! this also sets pointers to RESOURCES_SITE_REP
		! call resources_set_all_pointers(RESOURCES_EXCITON_REP)

		if (nr_blocks > 0) then

		!*************************************************************
		!  Calculation of the populations
		!*************************************************************
		if (resources_output_contains("populations")) then

			! output of populations

			call create_populations()
			have_populations = .true.


		end if

		!*************************************************************
		!  Calculation of coherences
		!*************************************************************
		if ((resources_output_contains(NOSE_RDM_B01) .or. &
		    resources_output_contains(NOSE_RDM_B01_ABS)).and. .not. have_coherences) then

			! output of coherences
			!if (external_data(1) == 1) then
			if (use_module_nakajima_zwanzig) then
				call read_coherences()
				have_coherences = .true.
			else
				call create_coherences()
				have_coherences = .true.
			end if

		end if

		!*************************************************************
		!  Calculation of linear polarization
		!*************************************************************
		if (resources_output_contains("polar_1")) then

			! output of the first order polarization

			!if (external_data(1) == 1 .and. .not. have_coherences) then
			if (use_module_nakajima_zwanzig .and. .not. have_coherences) then
				call read_coherences()
				have_coherences = .true.
			end if
		   	call create_polar_1()
		   	have_polar_1 = .true.

		end if

		!*************************************************************
		! Calculation of linear absorption spectrum
		!*************************************************************
		if (resources_output_contains("spect_abs")) then

			! output linear absorption spectrum

			!if (external_data(1) == 1) then
			if (use_module_nakajima_zwanzig) then
				if (.not.have_polar_1) then
					if (.not. have_coherences) then
						call read_coherences()
						have_coherences = .true.
					end if
					call create_polar_1()
				end if
			else
				if (.not.have_polar_1) then
					call create_polar_1()
				end if
			end if


		end if

		!*************************************************************
		!  Calculation of populations and non-optical coherences (NOC)
		!*************************************************************
		if ((resources_output_contains(NOSE_RDM_B11) .or. &
		resources_output_contains(NOSE_RDM_B11_ABS)).and. .not. have_pop_coh) then

			! output of populations and coherences

			!if (external_data(2) == 1) then
			if (use_module_nakajima_zwanzig) then
				call read_pop_coh()
				have_pop_coh = .true.
			else
				call create_pop_coh()
				call write_time_evolutions('E',.false.,.false.)
				have_pop_coh = .true.
			end if

		end if

		!*************************************************************
		!  Calculation of 2D spectrum
		!*************************************************************
		if (resources_output_contains("2d_ftpe")) then

			! output of 2D spectrum
			!if (external_data(3) == 1) then
			if (use_module_nakajima_zwanzig) then
		   		if (.not. have_coherences) then
					call read_coherences()
					have_coherences = .true.
				end if
				if (.not. have_pop_coh) then
					call read_pop_coh()
					have_pop_coh = .true.
				end if
				if (.not. have_1to2ex_coh) then
					call read_fe_coh()
					have_1to2ex_coh = .true.
				end if
			else
				if (.not. have_coherences) then
					call create_coherences()
					have_coherences = .true.
				end if
				if (.not. have_pop_coh) then
					call create_pop_coh()
					have_pop_coh = .true.
				end if
				if (.not. have_1to2ex_coh) then
					call create_fe_coh()
					have_1to2ex_coh = .true.
				end if
			end if



			call create_2D_spec()
			have_2d_ftpe = .true.

		end if

		!*************************************************************
		!  Calculation of 1 to 2 exciton coherences
		!*************************************************************
		if ((resources_output_contains(NOSE_RDM_B12) .or. &
		     resources_output_contains(NOSE_RDM_B12_ABS)).and. .not. have_1to2ex_coh) then

			! output of 1 to 2 exciton coherences

			!if (external_data(3) == 1) then
			if (use_module_nakajima_zwanzig) then
				call read_fe_coh()
				have_1to2ex_coh = .true.

			else
				call create_fe_coh()
				have_1to2ex_coh = .true.
			end if

		end if

		!*************************************************************
		!  Calculation of thermal light
		!*************************************************************
		if (resources_output_contains('thermal_response_function')) then

			if (.not. have_coherences) then

				if (use_module_nakajima_zwanzig) then
					call read_coherences()
					have_coherences = .true.
				else
					call create_coherences()
					have_coherences = .true.
				end if

			end if

			if (.not. have_pop_coh) then

				if (use_module_nakajima_zwanzig) then
					call read_pop_coh()
					have_pop_coh = .true.
				else
					call create_pop_coh()
					have_pop_coh = .true.
				end if

			end if

			!!!!!!!!!!!!!!!!!ACTIVE CODE WILL COME HERE!!!!!!!!!!!!!!!!!!!!!
			write(*,*) 'INITIALIZING DEBUG'
			call debug()

		end if

		end if ! if (nr_blocks > 0) then

	end subroutine do_qme_work


!-----------------------------------------------------------------------------------------------------------


	!
	! Calculates evolution of populations
	!
	subroutine create_populations()

		call weak_excitons_pop()

	end subroutine create_populations


	!
	! Collect populations from all processes
	!
	subroutine collect_populations

	end subroutine collect_populations


	!
	! Calculates evolution of the optical coherences
	!
	subroutine create_coherences()

		call weak_excitons_coh()

	end subroutine create_coherences


	!
	! Collect coherences from all processes
	!
	subroutine collect_coherences

	end subroutine collect_coherences


	!
	! Calculates evolution of NOC + coresponding populations
	!
	subroutine create_pop_coh()

		call weak_excitons_pop_coh()

	end subroutine create_pop_coh


	!
	! Calculates evolution of 1-to-2 exciton state coherences
	!
	subroutine create_fe_coh()

		call weak_excitons_fe_coh()

	end subroutine create_fe_coh


	!
	! Calculation of the first order polarization
	!
	subroutine create_polar_1()
		integer(i4b) :: i,j,k,b, mi,ni
		integer :: kk, kb
		! this will be a local pointer to the global transition dipole moments
		real(dp), dimension(:), pointer :: dd, dx,dy,dz
		! orientation averaging factors
		real(dp), dimension(:,:), allocatable :: as_orfact
		real(dp) :: vdni, vdmi

		if (.not.allocated(qme_polar_1)) then
			allocate(qme_polar_1(1:Nt(1)))
			qme_polar_1 = 0.0_dp
		end if

		if (.not.have_coherences) then

			call create_coherences()
			have_coherences = .true.

		end if

		call resources_rewind_blocks()

		! loop over blocks
		kk = 0
		kb = 1
		do

			if (use_twoexcitons) then

				dx => current_e_block%dx
		   		dy => current_e_block%dy
				dz => current_e_block%dz

			else

				dx => current_s_block%dx
				dy => current_s_block%dy
				dz => current_s_block%dz

			end if


			! this just temporary, it will be put somewhere else
			allocate(as_orfact(N1,N1))

			do ni = 1,N1
				do mi = 1,N1
					vdmi = sqrt(dx(mi)*dx(mi) + dy(mi)*dy(mi) + dz(mi)*dz(mi))
					vdni = sqrt(dx(ni)*dx(ni) + dy(ni)*dy(ni) + dz(ni)*dz(ni))
					if ((vdmi==0.0_dp).or.(vdni==0.0_dp)) then
					  as_orfact(ni,mi) = 0.0_dp
					else
					  as_orfact(ni,mi) = (1.0_dp/3.0_dp)* &
						(dx(ni)*dx(mi) + dy(ni)*dy(mi) + dz(ni)*dz(mi))/ &
						(vdni*vdmi)
					end if
				end do
			end do


!			allocate(dd(N1))
! pointing a pointer instead of allocation
			if (use_twoexcitons) then
				dd => current_e_block%dd
			else
				dd => current_s_block%dd
			end if
!			do i = 1, N1
!			   dd(i) = e(1)*dx(i) + e(2)*dy(i) + e(3)*dz(i)
!			end do

			! 4-mer test --- --- ---
			!write(*,*) '---> dd=', dd
			! 4-mer test --- --- ---

			! loop over time
			do i = 1, Nt(1)
				! loop over excited states
				do k = 1, N1
! this will work only in secular approximation, because only then all coherences
! start and end as the same ones. Otherwise they can transfer and we need to have
! and information about that. The information is given by evolution operator.
	 !			   qme_polar_1(i) = qme_polar_1(i) + (dd(k)**2)*as_orfact(k,k)*gcohs%C(i,k+kk)
!
! This how the full implementation should look like
! we need to implement the as_orfacts and Ueg
!
					do j = 1, N1
						qme_polar_1(i) = qme_polar_1(i) + (dd(k)*dd(j))*as_orfact(k,j)*evops(kb,kb)%Ueg(k,1,j,1,i) !gcohs%C(i,k+kk)

						! 4-mer test --- --- ---
						!write(*,*) '---> +++=', (dd(k)*dd(j))*as_orfact(k,j)*evops(kb,kb)%Ueg(k,1,j,1,i)
						!write(*,*) '---> +++=dddd  =', (dd(k)*dd(j))
						!write(*,*) '---> +++=orfact=', as_orfact(k,j)
						!write(*,*) '---> +++=evops =', evops(kb,kb)%Ueg(k,1,j,1,i)
						!write(*,*) '---> +++=kb	=', kb,' ',j,' ',k, ' ', i
						! 4-mer test --- --- ---

						!print *, "---> ", qme_polar_1(i), dd(k), dd(j), as_orfact(k,j), evops(kb,kb)%Ueg(k,1,j,1,i)
					end do
				end do
			end do

! no deallocation
!			deallocate(dd)
			deallocate(as_orfact)

			if (.not.resources_have_next_block()) exit

			kk = kk + N1
			kb = kb + 1
			call resources_next_block()


		end do

	end subroutine create_polar_1


	!
	! Collect polarizations from all processes
	!
	subroutine collect_polar_1()
	   integer :: N
		character(len=256) :: buff

		N = N_realizations_local
		call parallel_average_complex_array(qme_polar_1,N)
		qme_polar_1_collected = .true.

		write(buff,'(a,i5,a)') "Polar_1  : Total of ",N," realization(s) received"
		call print_log_message(trim(buff),5)

	end subroutine collect_polar_1


	!
	!  Calculates 2D spectrum
	!
	subroutine create_2D_spec()


		!if (parallel_id == 0) then
			if (.not.have_coherences) then
				call create_coherences()
				have_coherences = .true.
			end if

			if (.not.have_pop_coh) then
				call create_pop_coh()
				have_pop_coh = .true.
			end if

			if (.not.have_1to2ex_coh) then
				call create_fe_coh()
				have_1to2ex_coh = .true.
			end if

			call init_twod()
			if (.not.allocated(qme_2d_ftpe)) then
				allocate(qme_2d_ftpe(size(SR0,1),size(SR0,2)))
				qme_2d_ftpe = 0.0_dp
			 end if

			qme_2d_ftpe = qme_2d_ftpe + SR0
			qme_2d_ftpe = qme_2d_ftpe + SNR0

			call clean_twod()

		!end if


	end subroutine create_2D_spec


	!
	! Collects 2D spectra
	!
	subroutine collect_2D_spec()

		integer :: N
		character(len=256) :: buff

		N = N_realizations_local

		!call parallel_average_complex_array(polar_1,N)
		!polar_1_collected = .true.

		!write(buff,'(a,i5,a)') "Polar_1  : Total of ",N," realization(s) received"
		!call print_log_message(trim(buff),5)

		!qme_2d_ftpe = qme_2d_ftpe/N_realizations_local
		call parallel_average_complex_array2(qme_2d_ftpe,N)

		write(buff,'(a,i5,a)') "2d_ftpe  : Total of ",N," realization(s) received"
		call print_log_message(trim(buff),5)

	end subroutine collect_2D_spec


	!
	! Cleans all allocated memory
	!
	subroutine clean_qme

		call clean_resources_qme()

	end subroutine clean_qme


!-----------------------------------------------------------------------------------------------------------


	!
	! Collecting data after simulation finished
	!
	subroutine collect_qme_data(err)
		integer, intent(out) :: err

		integer :: k,i, fs
		real(dp) :: oma

		character(len = 20) 	:: buffer
		character(len = 300)	:: name
		character(len = 80) 	:: caption_line ! with intention of proceding output data with Origin

		err = 0

		!*************************************************************
		!  Outputting populations
		!*************************************************************
		if (resources_output_contains("populations")) then

			! output of the populations

			call collect_populations()

			if (parallel_id == 0) then
				open(unit=11,file=trim(file_join(out_dir,"populations.dat")))

			  !-------------------- Organizes the caption line for Origin
			  caption_line = 'time'
			  do i = 1, pops%N
				write(buffer,'(i1)') i
				buffer = trim('Population')//trim(buffer)
				caption_line = trim(caption_line)//' '//trim(buffer)
				print *, buffer
			  end do
			  write(11,'(a, /)') caption_line
			  !--------------------


				write(buffer,'(i3)') pops%N+1
				buffer = '('//trim(buffer)//'f12.8)' ! forming a descriptor for outputting row by row

				do i = 1, Nt(1)

				  write(11,buffer) (i-1)*gt(1)*dt, pops%P(i,1:pops%N)

				end do


				close(unit=11)
			end if

		end if

		!*************************************************************
		!  Outputting coherences
		!*************************************************************
		if (resources_output_contains("coherences")) then

			! output of the coherences

			call collect_coherences()

			if (parallel_id == 0) then
				open(unit=11,file=trim(file_join(out_dir,"coherences.dat")))

!				 do k = 1, gcohs%N
!				 do i = 1, Nt(1)
!
!					write(11,*) (i-1)*gt(1)*dt, real(gcohs%C(i,k)), aimag(gcohs%C(i,k))
!
!				 end do
!				 end do

				write(buffer,'(i3)') gcohs%N*2+1
				buffer = '('//trim(buffer)//'f12.8)' ! forming a descriptor for outputting row by row

				do i = 1, Nt(1)

				  write(11,buffer) (i-1)*gt(1)*dt, real(gcohs%C(i,1:gcohs%N)), aimag(gcohs%C(i,1:gcohs%N))

				end do

				close(unit=11)
			end if

		end if

		!*************************************************************
		!  Outputting linear polarization
		!*************************************************************
		if (resources_output_contains("polar_1")) then

			! output of the polarization

			call collect_coherences()

			if (parallel_id == 0) then
				open(unit=11,file=trim(file_join(out_dir,"polar_1.dat")))

				 do i = 1, Nt(1)

					write(11,*) (i-1)*gt(1)*dt, real(qme_polar_1(i)), aimag(qme_polar_1(i))

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

			!oma = maxval(qme_spect_abs)
			!qme_spect_abs = qme_spect_abs/oma

			if (parallel_id == 0) then
				open(unit=11,file=trim(file_join(out_dir,"spect_abs.dat")))

				do i = 1, NFFT
					 oma = ((-NFFT/2 + i)*dom + rwa)*Energy_internal_to_cm

					write(11,*) oma, qme_spect_abs(i)

				end do


				close(unit=11)
			end if

		end if

		!*************************************************************
		!  Outputting NOC's + populations
		!*************************************************************
		if (resources_output_contains(NOSE_RDM_B11)) then


			if (parallel_id == 0) then
				open(unit=11,file=trim(file_join(out_dir,"pop_n_o_coh.dat")))

				write(buffer,'(i3)') 2*((rcohs%N)**2) + 1
				buffer = '('//trim(buffer)//'f16.8)'

				do i = 1, Nt(1)

				  write(11,buffer) (i-1)*gt(1)*dt, evops(1,1)%Uee(1:rcohs%N,1:rcohs%N,1,1,i)
				  !rcohs%RC(i,1:rcohs%N,1:rcohs%N)

				end do

				close(unit=11)

			end if

		end if

		!*************************************************************
		!  Outputting 2D spectra
		!*************************************************************
		if (resources_output_contains("2d_ftpe")) then

			call collect_2D_spec()

			if (parallel_id == 0) then

				fs = gt(2)*dt
				if (fs < 10) then
					write(buffer,'(i1)') fs
				elseif (fs < 100) then
					write(buffer,'(i2)') fs
				elseif (fs < 1000) then
					write(buffer,'(i3)') fs
				else
					write(buffer,'(i4)') fs
				endif


				name = trim('twod_re_tot_T=')//trim(buffer)//trim('fs.dat')
				call save_2D_tot(real(qme_2d_ftpe),trim(file_join(out_dir,name)))
				name = trim('twod_im_tot_T=')//trim(buffer)//trim('fs.dat')
				call save_2D_tot(aimag(qme_2d_ftpe),trim(file_join(out_dir,name)))

                ! temporarily outputing limits from here
                call save_2D_limits(dom,NFFT,trim(file_join(out_dir,"limits.dat")))

            end if

		end if

		!*************************************************************
		!  Outputting 1- to 2- exciton coherences
		!*************************************************************
		if (resources_output_contains(NOSE_RDM_B12)) then

			if (parallel_id == 0) then
				open(unit=11,file=trim(file_join(out_dir,"1to2_ex_coh.dat")))

				write(buffer,'(i3)') 2*2 + 1
				buffer = '('//trim(buffer)//'f16.8)'

				do i = 1, Nt(1)

				  write(11,buffer) (i-1)*gt(1)*dt, evops(1,1)%Ufe(1,1:2,1,1,i)

				end do

				close(unit=11)

			end if

		end if


	end subroutine collect_qme_data


!-----------------------------------------------------------------------------------------------------------


	!
	! Calculate linear absorption spectrum from polarization
	!
	subroutine create_spect_abs()

		complex(dpc), dimension(:),   allocatable	:: sig
		complex(dpc), dimension(:,:), allocatable	:: dat
		real(dp), dimension(:), allocatable :: rs

		integer :: on, ch, i
		real(dp) :: oma, rr

		padfac = 4
		NFFT = (2**padfac)*NFFT_basic

		allocate(sig(NFFT))
		allocate(dat(1,NFFT))
		allocate(rs(NFFT))
		if(.not. allocated(qme_spect_abs)) then
			allocate(qme_spect_abs(NFFT))
		end if
		rs			= 0.0_dp
		qme_spect_abs = 0.0_dp

		dom = 2.0_dp*PI_D/(NFFT*gt(1)*dt)

		if (.not.qme_polar_1_collected) then
			call collect_polar_1()
		end if


		sig = 0.0_dp
		sig(1:Nt(1)) = qme_polar_1(1:Nt(1))

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

!				if (i == 0) then
!					rs(i) = abs(oma)*real(sig(NFFT))
!				else
					rs(i) = abs(oma)*real(sig(i))
!				end if


		end do

		rr = 1 !maxval(rs(:))

		qme_spect_abs = rs !/rr
		qme_spect_abs = qme_spect_abs * (gt(1)*dt)

		deallocate(sig)
		deallocate(dat)
		deallocate(rs)

	end subroutine create_spect_abs

	subroutine read_coherences()
		integer(i4b)	  :: N2
		call print_log_message("Reading coherences",5)
		N2 = N1*(N1-1)/2

		!call debug()

		!call read_data(1,N1,1, N1, 'O', submethod1, Nt(1)*gt(1))
		call fill_evolution_superoperator_nakajima_zwanzig('O',submethod1)
		if (resources_output_contains(NOSE_RDM_B01)) then
			call write_time_evolutions('O',.false.,.false.)
		endif
		if (resources_output_contains(NOSE_RDM_B01_ABS)) then
			call write_time_evolutions('O',.true.,.false.)
		endif
		if (resources_output_contains(NOSE_RDM_B01_CONJG)) then
			call write_time_evolutions('O',.false.,.true.)
		endif
		call write_evolution_operators('O')
		call write_redfield_tensor('O')


	end subroutine read_coherences

	subroutine read_pop_coh()
		integer(i4b) :: k, N2
		call print_log_message("Reading populations-coherences",5)

		!call read_data(1,N1*N1,1, N1*N1, 'E', submethod2, Nt(1)*gt(1))
		call fill_evolution_superoperator_nakajima_zwanzig('E',submethod2)
		if (resources_output_contains(NOSE_RDM_B11)) then
			call write_time_evolutions('E',.false.,.false.)
		endif
		if (resources_output_contains(NOSE_RDM_B11_ABS)) then
			call write_time_evolutions('E',.true.,.false.)
		endif
		if (resources_output_contains(NOSE_RDM_B11_CONJG)) then
			call write_time_evolutions('E',.false.,.true.)
		endif

	end subroutine read_pop_coh

	subroutine read_fe_coh()
		integer(i4b) :: k, N2
		k = 1
		call print_log_message("Reading 1-2 coherences",5)

		!call read_data(1,N1*N1*(N1-1)/2,1, N1*N1*(N1-1)/2, '2', submethod1, Nt(1)*gt(1))
		call fill_evolution_superoperator_nakajima_zwanzig('2',submethod1)
		if (resources_output_contains(NOSE_RDM_B12)) then
			call write_time_evolutions('2',.false.,.false.)
		endif
		if (resources_output_contains(NOSE_RDM_B12_ABS)) then
			call write_time_evolutions('2',.true.,.false.)
		endif
		if (resources_output_contains(NOSE_RDM_B12_CONJG)) then
			call write_time_evolutions('2',.false.,.true.)
		endif

	end subroutine read_fe_coh


	!*************************************************************
	!  Writing out density matrix evolution after delta-excitation
	!*************************************************************

	subroutine write_time_evolutions(type,absolute_value,conjugate_coherences)
		! writes time evolution of the density matrix elements according to the type,
		! type = 'E', 'O', '2'. If type == '2', arrays are allocated as 0:N1+N1*(N-1)/2,
		! otherwise 0:N1

      	character, intent(in)	:: type
      	logical	, intent(in)	:: absolute_value, conjugate_coherences

		real(dp), dimension(:,:), allocatable 			:: dd
      	real(dp) 											:: s, x, magnitude
      	complex(dpc)										:: s_c
      	integer 											:: a,b,i,j,k,l,Ublock
      	complex(dpc), dimension(N1_from_type(type),N2_from_type(type)) 	:: rho0, eig1, eig2
      	complex(dpc), dimension(:,:,:), allocatable 						:: rr
     	complex(dpc), dimension(:,:,:,:,:), pointer		:: actual_U
      	character(len=10)	:: cha,chb
      	character(len=64)	:: name

      	if(.not. (type == 'E' .or. type == 'O' .or. type == '2')) then
 			call print_error_message(-1,  'type error in write_time_evolutions')
 			stop
      	end if

!		open(unit=13,file='/home/olsij4am/prace/nose-analytical-coherence-markov.dat')
!		open(unit=14,file='/home/olsij4am/prace/nose-analytical-coherence-tau.dat')

!		write(*,*) 'writing debug functions - coherence from cummulant expansion'
!
!		do i=1, Nt(1)
!			if((i-1)*dt*gt(1)+tau_of_projector >= gt(1)*dt*Nt(1)) then
!				exit
!			end if
!
!			write(13,*) (i-1)*dt*gt(1), &
!				abs(current_e_block%dd(1)**2)*abs(exp(-conjg(goft_exciton(2,2,(i-1)*dt*gt(1))) -goft_exciton(1,1,(i-1)*dt*gt(1)) -conjg(goft_exciton(2,2,tau_of_projector))  )), &
!				abs(current_e_block%dd(1)**2)*real(exp(-conjg(goft_exciton(2,2,(i-1)*dt*gt(1))) -goft_exciton(1,1,(i-1)*dt*gt(1)) -conjg(goft_exciton(2,2,tau_of_projector))  )), &
!				abs(current_e_block%dd(1)**2)*aimag(exp(-conjg(goft_exciton(2,2,(i-1)*dt*gt(1))) -goft_exciton(1,1,(i-1)*dt*gt(1)) -conjg(goft_exciton(2,2,tau_of_projector))  ))
!!				real(goft_exciton(1,1,(i-1)*dt*gt(1))),real(goft_exciton(2,2,(i-1)*dt*gt(1)))
!			write(14,*) (i-1)*dt*gt(1), &
!				abs(current_e_block%dd(1)**2)*abs(exp(-conjg(goft_exciton(2,2,(i-1)*dt*gt(1)+tau_of_projector)) - goft_exciton(1,1,(i-1)*dt*gt(1))  )), &
!				abs(current_e_block%dd(1)**2)*real(exp(-conjg(goft_exciton(2,2,(i-1)*dt*gt(1)+tau_of_projector)) - goft_exciton(1,1,(i-1)*dt*gt(1))  )), &
!				abs(current_e_block%dd(1)**2)*aimag(exp(-conjg(goft_exciton(2,2,(i-1)*dt*gt(1)+tau_of_projector)) - goft_exciton(1,1,(i-1)*dt*gt(1))  ))
!
!		end do
!
!		close(13)
!		close(14)


		Ublock = 1

		allocate(rr(N1_from_type(type),N2_from_type(type),Nt(1)))
      	rr = 0.0d0

		! We set indices range according to block we evaluate. Because rho0 is
		! whole density matrix, while evolution operators are only from particular
		! block, offset is set between these indices.
		if (type == '2') then
			actual_U => evops(Ublock,Ublock)%Ufe
		else if (type == 'E') then
			actual_U => evops(Ublock,Ublock)%Uee
		else if (type == 'O') then
			actual_U => evops(Ublock,Ublock)%Ueg
		end if

		call calculate_dipoled_rho(type,rho0)

      	do a = 1, N1_from_type(type)
      	do b = 1, N2_from_type(type)

      	write(cha,'(i1)') a
      	write(chb,'(i1)') b

		if (type == 'E') then
			name = 'dens_E_'//trim(cha)//'_'//trim(chb)//'.dat'
      	else if (type == 'O') then
      		name = 'dens_O_'//trim(cha)//'_'//trim(chb)//'.dat'
      	else if (type == '2') then
      		name = 'dens_2_'//trim(cha)//'_'//trim(chb)//'.dat'
      	end if

      	if(absolute_value .eqv. .true.) then
      		name = 'abs_' // name
      	endif

      	if(conjugate_coherences .eqv. .true.) then
      		name = 'conjg_' // name
      	endif

		open(unit=11,file=trim(file_join(out_dir,trim(name))))

      	do i = 1, Nt(1)

	  		do k = 1, N1_from_type(type)
			do l = 1, N2_from_type(type)

				if(conjugate_coherences .eqv. .true.) then
					!!! averaging coherences
     				rr(a,b,i) = rr(a,b,i) + (actual_U(a,b,k,l,i)*rho0(k,l) 		  &
     									  + conjg(actual_U(b,a,k,l,i)*rho0(k,l) ) &
     									   )/2

				else
     				rr(a,b,i) = rr(a,b,i) + actual_U(a,b,k,l,i)*rho0(k,l)
     			end if

       		end do
     		end do

			if (absolute_value .eqv. .true.) then
				write(11,*) (i-1)*dt*gt(1), abs(rr(a,b,i))
			else
				write(11,*) (i-1)*dt*gt(1), real(rr(a,b,i)), aimag(rr(a,b,i))
          	endif


      	end do

      	close(unit=11)

      	end do
      	end do

      	if (type == 'E') then
      		name = 'dens_trace.dat'
      		open(11,file=trim(file_join(out_dir,trim(name))))

      		do i = 1, Nt(1)
         		s_c = 0.0d0
         		do a = 1, N1
            		s_c = s_c + rr(a,a,i)
         		end do

		 		write(11,*) (i-1)*dt*gt(1), real(s_c), aimag(s_c)
    	  	end do

    	  	close(11)

    	  	do a = 1, N1_from_type(type)
    	  		write(cha,'(i1)') a
	    	  	name = 'dens_E_eigenvalue_'//trim(cha)//'.dat'
    	  		open(11,file=trim(file_join(out_dir,trim(name))))

      			do i = 1, Nt(1)
					call spec(rr(:,:,i),eig1,eig2)

			 		write(11,*) (i-1)*dt*gt(1), real(eig2(a,a)), aimag(eig2(a,a))
    	  		end do

    	  		close(11)
    	  	end do
		end if


      	deallocate(rr)

	end subroutine write_time_evolutions

!	Former version of previous output function. Will be deleted after proper testing of
!		previdous one. In case of problems with output, try to replace new function by
!		this one...
!
!	subroutine write_time_evolutions(type,absolute_value)
!		! writes time evolution of the density matrix elements according to the type,
!		! type = 'E', 'O', '2'. If type == '2', arrays are allocated as 0:N1+N1*(N-1)/2,
!		! otherwise 0:N1
!
!		double complex, dimension(:), allocatable 		:: rho
!		integer 											:: N, N2
!		integer												:: FileMaxIndex1, FileMaxIndex2
!		integer												:: BlockOffset1, BlockOffset2
!		double precision, dimension(:,:), allocatable 	:: dd
!     	double precision 									:: s, x, magnitude
!    	integer 											:: a,b,i,j,k,l,Ublock
!      	double complex, dimension(:,:,:), allocatable 	:: rr
!     	double complex, dimension(:,:), allocatable 	:: rho0
!      	complex(dpc), dimension(:,:,:,:,:), pointer		:: actual_U
!     	character(len=10)	:: cha,chb
!    	character(len=64)	:: name
!   	character			:: type
!      	logical				:: absolute_value
!
!     	!write(*,*) typ, N1, Nt(1)
!    	if(.not. (type == 'E' .or. type == 'O' .or. type == '2')) then
!   		return
!     	end if
!
!		Ublock = 1
!     	N2 = N1*(N1-1)/2
!		N = N1 + N1*(N1-1)/2
!
!     	allocate(rho(Nt(1)),rr(N,N,Nt(1)),rho0(0:N,0:N))
!     	allocate(dd(0:N,0:N))
!
!      	rho0 = 0.0d0
!    	rr = 0.0d0
!     	dd = 0.0d0
!
!     	rho0(0,0) = 1.0d0
!
!     	do k = 1, N1
!      		!write(*,*) 'current_e_block%dd(k)=',current_e_block%dd(k)
!  			dd(0,k) = current_e_block%dd(k)
!  			dd(k,0) = current_e_block%dd(k)
!     	end do
!
!		! we add "2-coherence" block of dd if 2-coherences are computed
!      	if (type == '2') then
!     		do k = 1, N1
!			do l = 1, N2
!  				dd(k,l+N1) = current_e_block%dd_2(l,k)
!  				dd(l+N1,k) = current_e_block%dd_2(l,k)
!  				!write(*,*) 'current_e_block%dd(k)=',current_e_block%dd(k)
!     		end do
!      		end do
!     	end if
!
!     	close(12)
!
!		!if(type == '2') then
!		!	write(*,*) 'dd  =',dd
!		!endif
!
!     	x = sum(dd(0,:)**2)
!
!		! normalizace na jednu desetinu? pridana odmocnina
!		magnitude = 1.0d0/10
!     	dd = dd/(sqrt(x)) * magnitude
!
!
!		if (type == 'E') then
!			! only this term is relevant for the Uee block, 2 in commutator cancels with
!			! the 2 from the Taylor expansion, other terms fo to the other blocks of U
!     		rho0 = 	matmul(dd,matmul(rho0,dd))
!      	else if (type == 'O') then
!     		rho0 = 	cmplx(0,1)*(matmul(dd,rho0)-matmul(rho0,dd))
!      	else if (type == '2') then
!     		rho0 = 	cmplx(0,1)*(matmul(dd,matmul(dd,matmul(rho0,dd)))- &
!      				matmul(matmul(dd,matmul(rho0,dd)),dd))
!     	end if
!
!		!if(type == '2') then
!		!	write(*,*) 'rho0=',rho0
!		!endif
!
!		! We set indices range according to block we evaluate. Because rho0 is
!		! whole density matrix, while evolution operators are only from particular
!		! block, offset is set between these indices.
!		if (type == '2') then
!			FileMaxIndex1 = N2
!			FileMaxIndex2 = N1
!			BlockOffset1  = N1
!			BlockOffset2  = 0
!			actual_U => evops(Ublock,Ublock)%Ufe
!		else if (type == 'E') then
!			FileMaxIndex1 = N1
!			FileMaxIndex2 = N1
!			BlockOffset1  = 0
!			BlockOffset2  = 0
!			actual_U => evops(Ublock,Ublock)%Uee
!		else if (type == 'O') then
!			FileMaxIndex1 = N1
!			FileMaxIndex2 = 1
!			BlockOffset1  = 0
!			BlockOffset2  = -1
!			actual_U => evops(Ublock,Ublock)%Ueg
!		end if
!
!      	do a = 1, FileMaxIndex1
!     	do b = 1, FileMaxIndex2
!
!      	write(cha,'(i1)') a
!     	write(chb,'(i1)') b
!
!		if (type == 'E') then
!			name = 'dens_E_'//trim(cha)//'_'//trim(chb)//'.dat'
!     	else if (type == 'O') then
!      		name = 'dens_O_'//trim(cha)//'_'//trim(chb)//'.dat'
!     	else if (type == '2') then
!      		name = 'dens_2_'//trim(cha)//'_'//trim(chb)//'.dat'
!     	end if
!
!     	if(absolute_value .eqv. .true.) then
!      		name = 'abs_' // name
!     	endif
!
!		open(unit=11,file=trim(file_join(out_dir,trim(name))))
!
!     	do i = 1, Nt(1)
!
!	  		do k = 1, FileMaxIndex1
!				do l = 1, FileMaxIndex2
!	     			rr(a,b,i) = rr(a,b,i) + actual_U(a,b,k,l,i)*rho0(BlockOffset1+k,BlockOffset2+l)
!	     			!write(*,*) 'U=', actual_U(a,b,k,l,i)
!	     			!write(*,*) 'rho0=', rho0(BlockOffset1+k,BlockOffset2+l), k, l
!	     			!write(*,*) 'rr=', rr(a,b,i)
!
!          		end do
!    		end do
!
!			if (absolute_value .eqv. .true.) then
!				write(11,*) (i-1)*dt*gt(1), abs(rr(a,b,i))
!			else
!          		write(11,*) (i-1)*dt*gt(1), real(rr(a,b,i)), aimag(rr(a,b,i))
!         	endif
!
!
!     	end do
!
!     	close(11)
!
!      	end do
!     	end do
!
!     	if (type == 'E') then
!      		name = 'dens_trace.dat'
!     		open(11,file=trim(file_join(out_dir,trim(name))))
!
!      		do i = 1, Nt(1)
!        		s = 0.0d0
!         		do a = 1, N
!	          		s = s + rr(a,a,i)
!        		end do
!
!		 		write(11,*) real(s)
!   	  	end do
!		end if
!
!     	close(11)
!
!
!     	deallocate(rho,rr,rho0)
!
!	end subroutine write_time_evolutions






	!*************************************************************
	!  Calculate rho times \mu specific for delta-pulse excit.
	!*************************************************************

	subroutine calculate_dipoled_rho(type_in,rho_out)

		character, intent(in)							:: type_in
		complex(dpc), dimension(:,:), intent(out) 	:: rho_out

		complex(dpc), dimension(N1_from_type('O'),N2_from_type('O')) :: rho_O
		complex(dpc), dimension(N1_from_type('E'),N2_from_type('E')) :: rho_E
		complex(dpc), dimension(N1_from_type('g'),N2_from_type('g')) :: rho_g

		integer :: i,j
		logical :: get_tau_phase

     	if(.not. (type_in == 'E' .or. type_in == 'O' .or. type_in == '2')) then
			call print_error_message(-1,  'type_in error in calculate_dipoled_rho')
			stop
    	end if

      	if(.not. (size(rho_out,1) == N1_from_type(type_in) .and. size(rho_out,2) == N2_from_type(type_in))) then
 			call print_error_message(-1,  'dimension error2 in calculate_dipoled_rho')
 			stop
      	end if



		rho_g(1,1) = 1
		get_tau_phase = .true.

		if(type_in == 'O') then
			call calculate_dipole_excitation('g','O',rho_g,rho_out,.true.)
		else if(type_in == 'E') then
			call calculate_dipole_excitation('g','O',rho_g,rho_O,.true.)

			if(get_tau_phase) then
				do i=1,N1_from_type('O')
				do j=1,N2_from_type('O')
					rho_O(i,j) = rho_O(i,j)*exp(cmplx(0,1,dpc)*tau_of_projector &
							*(iblocks(1,1)%eblock%en(j) - rwa) 					&
							- goft_exciton(j,j,tau_of_projector))
				end do
				end do
			end if

			call calculate_dipole_excitation('O','E',rho_O,rho_out,.false.)

			if(submethod2 /= 'u' .and. submethod2 /= 'U' .and. .not. tau_projector_normalization_for_others) then
				rho_out = 2*rho_out
			else
				! normalization to trace 1 for tau = 0
				call calculate_dipole_excitation('g','O',rho_g,rho_O,.true.)
				call calculate_dipole_excitation('O','E',rho_O,rho_E,.false.)

				rho_out = rho_out / trace(rho_E)
			end if

!			if(get_tau_phase) then
!				do i=1,N1_from_type('E')
!				do j=1,N2_from_type('E')
!					rho_out(i,j) = rho_out(i,j)*exp(cmplx(0,1,dpc)*tau_of_projector &
!							*(iblocks(1,1)%eblock%en(j) - rwa) 					&
!							- goft_exciton(j,j,tau_of_projector))
!				end do
!				end do
!			end if

		else if(type_in == '2') then
			call calculate_dipole_excitation('g','O',rho_g,rho_O,.true.)
			call calculate_dipole_excitation('O','E',rho_O,rho_out,.true.)
		end if

	end subroutine calculate_dipoled_rho


	!*************************************************************
	!  Calculate rho times \mu
	!*************************************************************

	subroutine calculate_dipole_excitation(type_in,type_out,rho_in,rho_out,mu_from_left)
		! writes time evolution of the density matrix elements according to the type,
		! type = 'E', 'O', '2'. If type == '2', arrays are allocated as 0:N1+N1*(N-1)/2,
		! otherwise 0:N1

      	character, intent(in)							:: type_out
      	character, intent(in)							:: type_in
    	complex(dpc), dimension(:,:), intent(in) 		:: rho_in
      	complex(dpc), dimension(:,:), intent(out) 	:: rho_out
      	logical, intent(in)							:: mu_from_left

		integer 											:: N, N2
!		integer												:: FileMaxIndex1, FileMaxIndex2
!		integer												:: BlockOffset1, BlockOffset2
		real(dp), dimension(:,:), allocatable 			:: dd
      	real(dp) 											:: s, x, magnitude
      	!complex(dpc)										:: s_c
      	integer 											:: a,b,i,j,k,l,Ublock
      	complex(dpc), dimension(:,:), allocatable 		:: rho0

      	if(.not. (type_out == 'E' .or. type_out == 'O' .or. type_out == '2'.or. type_out == 'g')) then
 			call print_error_message(-1,  'type error in calculate_dipole_excitation')
 			stop
      	end if

     	if(.not. (type_in == 'E' .or. type_in == 'O' .or. type_in == '2' .or. type_in == 'F' .or. type_in == 'g')) then
			call print_error_message(-1,  'type_in error in calculate_dipole_excitation')
			stop
    	end if

      	if(.not. (size(rho_out,1) == N1_from_type(type_out) .and. size(rho_out,2) == N2_from_type(type_out))) then
 			call print_error_message(-1,  'dimension error1 in calculate_dipole_excitation')
 			stop
      	end if

      	if(.not. (size(rho_in,1) == N1_from_type(type_in) .and. size(rho_in,2) == N2_from_type(type_in))) then
 			call print_error_message(-1,  'dimension error2 in calculate_dipole_excitation')
 			stop
      	end if

		Ublock = 1
      	N2 = N1*(N1-1)/2
		N = N1 + N1*(N1-1)/2

      	allocate(rho0(0:N,0:N))
      	allocate(dd(0:N,0:N))

      	rho0 = 0.0d0
      	rho_out = 0.0d0
      	dd = 0.0d0

		! put rho_in inside bigger matrix
		if(type_in == 'g') then
      		rho0(0,0) = rho_in(1,1)
		else if(type_in == 'O') then
			rho0(1:N1,0:0) = rho_in
		else if(type_in == 'E') then
			rho0(1:N1,1:N1) = rho_in
		else if(type_in == '2') then
			rho0((N1+1):N,1:N1) = rho_in
		else if(type_in == 'F') then
		end if

		! create bigger matrix of dipole moments
		if (use_module_nakajima_zwanzig) then
      		do k = 1, N1
   				dd(0,k) = current_e_block%dd(k)
   				dd(k,0) = current_e_block%dd(k)
      		end do
      	else
      		do k = 1, N1
   				dd(0,k) = current_s_block%dd(k)
   				dd(k,0) = current_s_block%dd(k)
      		end do
      	end if

		! we add "2-coherence" block of dd if 2-coherences are computed
      	if (type_out == '2') then
      		do k = 1, N1
			do l = 1, N2
   				dd(k,l+N1) = current_e_block%dd_2(l,k)
   				dd(l+N1,k) = current_e_block%dd_2(l,k)
   				!write(*,*) 'current_e_block%dd(k)=',current_e_block%dd(k)
      		end do
      		end do
      	end if

 !     	x = sum(dd(0,:)**2)

		! normalizace na jednu desetinu? pridana odmocnina
!		magnitude = 1.0d0/10
!      	dd = dd/(sqrt(x)) * magnitude

		if (mu_from_left) then
      		rho0 = 	cmplx(0,1,dpc)*matmul(dd,rho0)
		else
			rho0 = 	-cmplx(0,1,dpc)*matmul(rho0,dd)
      	end if

		! We set indices range according to block we evaluate. Because rho0 is
		! whole density matrix, while evolution operators are only from particular
		! block, offset is set between these indices.
		if(type_out == 'g') then
      		rho_out(1,1) = rho0(0,0)
		else if(type_out == 'O') then
			rho_out = rho0(1:N1,0:0)
		else if(type_out == 'E') then
			rho_out = rho0(1:N1,1:N1)
		else if(type_out == '2') then
			rho_out = rho0((N1+1):N,1:N1)
		end if

      	deallocate(rho0,dd)

	end subroutine calculate_dipole_excitation

	function single_photon_I_laser(t1,t2,t0,Delta) result(ret)
		real(dp), intent(in)	:: t1, t2, t0, Delta
		real(dp)				:: ret

		ret = exp(-(t1-t0)*(t1-t0)/Delta/Delta-(t2-t0)*(t2-t0)/Delta/Delta)
	end function single_photon_I_laser

	function single_photon_I_thermal(t1,t2,tau,w) result(ret)
		real(dp), intent(in)	:: t1, t2, tau, w
		real(dp)				:: ret

		ret = exp(-abs(t1-t2)/tau + cmplx(0,w,dpc)*(t1-t2))
	end function single_photon_I_thermal


	subroutine integrate_pumped_rho_laser(type_pump,rho_pump,t0,Delta,T_perioda)
		! integrates excitation of rho by light of type type_I
		! 't' -- thermal, 'l' -- periodic laser pulses

      	character, intent(in)							:: type_pump
      	complex(dpc), dimension(:,:,:), intent(out) 	:: rho_pump
      	real(dp), intent(in)							:: Delta, t0, T_perioda

      	real(dp)	:: time, ttime, tttime
      	integer		:: i0,i,ii,iii, k,l,r,s,u,v

		complex(dpc), dimension(N1_from_type('g'),N2_from_type('g'))	:: rho_g
      	complex(dpc), dimension(N1_from_type('O'),N2_from_type('O'))	:: rho_O, rho_O2
      	complex(dpc), dimension(N1_from_type('E'),N2_from_type('E'))	:: rho_E, tmp1, tmp2
      	complex(dpc), dimension(size(rho_pump,1),size(rho_pump,2),size(rho_pump,3)) &
      					:: drho_pump


      	if(.not. (type_pump == 'E' .or. type_pump == 'O' .or. type_pump == '2'.or. type_pump == 'g')) then
 			call print_error_message(-1,  'type error in integrate_pumped_rho')
 			stop
      	end if

      	if(.not. (size(rho_pump,1) == N1_from_type(type_pump) .and. size(rho_pump,2) == N2_from_type(type_pump))) then
 			call print_error_message(-1,  'dimension error1 in integrate_pumped_rho')
 			stop
      	end if


		rho_pump = 0.0_dp
		rho_g(1,1) = 1.0_dp
		call calculate_dipole_excitation('g','O',rho_g,rho_O,.true.)

!				call calculate_dipole_excitation('g','O',rho_g,rho_O,.true.)
!				call calculate_dipole_excitation('O','E',rho_O,rho_E,.true.)
!				write(*,*) 'TT: ',rho_E
!				call calculate_dipole_excitation('g','O',rho_g,rho_O,.true.)
!				call calculate_dipole_excitation('O','E',rho_O,rho_E,.false.)
!				write(*,*) 'TF: ',rho_E
!				call calculate_dipole_excitation('g','O',rho_g,rho_O,.false.)
!				call calculate_dipole_excitation('O','E',rho_O,rho_E,.true.)
!				write(*,*) 'FT: ',rho_E
!				call calculate_dipole_excitation('g','O',rho_g,rho_O,.false.)
!				call calculate_dipole_excitation('O','E',rho_O,rho_E,.false.)
!				write(*,*) 'FF: ',rho_E
!				stop

		i0 = 1
      	do i = i0, 150 ! Nt(1)
     		write(*,*) 'INTEGRUJI CAS T = ',(i-1)*dt*gt(1),' fs'
      	do ii = i0, i
      	do iii = i0, ii

      		ttime  = (ii-1)*dt*gt(1)
      		tttime = (iii-1)*dt*gt(1)

      		rho_O2 = 0.0_dp

			do k=1,N1_from_type('O')
			do l=1,N2_from_type('O')

			do r=1,N1_from_type('O')
			do s=1,N2_from_type('O')

				rho_O2(k,l) = rho_O2(k,l) + evops(1,1)%Ueg(k,l,r,s, ii-iii+1)*rho_O(r,s)

			end do
			end do

			end do
			end do

			call calculate_dipole_excitation('O','E',rho_O2,rho_E,.false.)

			do k=1,N1_from_type('E')
			do l=1,N2_from_type('E')

			do r=1,N1_from_type('E')
			do s=1,N2_from_type('E')

			! HLOUPA IMPLEMENTACE PERIODY -- TAKTO TRVA O(3) -- JE TREBA SPOCITAT
			! DERIVACI RHO A PAK JI ZPERIODIZOVAT...
			do u = 0, INT(2*max(ttime,tttime)/T_perioda)

				rho_pump(k,l,i) = rho_pump(k,l,i) 	 								&
					+ evops(1,1)%Uee(k,l,r,s,i-ii+1)*rho_E(r,s)						&
					* single_photon_I_laser(ttime,tttime,t0+T_perioda*u,Delta)		&
					* gt(1)*dt * gt(1)*dt

			write(*,*) single_photon_I_laser(ttime,tttime,t0+T_perioda*u,Delta)

			end do

			end do
			end do

			end do
			end do

      	end do
      	end do

      	rho_pump(:,:,i) = rho_pump(:,:,i) + transpose(conjg(rho_pump(:,:,i)))

      	end do

	end subroutine integrate_pumped_rho_laser

	subroutine integrate_pumped_rho_thermal(type_pump,rho_pump,tau,w)
		! integrates excitation of rho by thermal light

      	character, intent(in)							:: type_pump
      	complex(dpc), dimension(:,:,:), intent(out) 	:: rho_pump
      	real(dp)										:: tau, w

      	real(dp)	:: time, ttime, tttime
      	integer		:: i0,i,jj,jjj, k,l,r,s

		complex(dpc), dimension(N1_from_type('g'),N2_from_type('g'))	:: rho_g
      	complex(dpc), dimension(N1_from_type('O'),N2_from_type('O'))	:: rho_O, rho_O2
      	complex(dpc), dimension(N1_from_type('E'),N2_from_type('E'))	:: rho_E, tmp1, tmp2
      	complex(dpc), dimension(size(rho_pump,1),size(rho_pump,2),size(rho_pump,3)) &
      					:: drho_pump


      	if(.not. (type_pump == 'E' .or. type_pump == 'O' .or. type_pump == '2'.or. type_pump == 'g')) then
 			call print_error_message(-1,  'type error in integrate_pumped_rho')
 			stop
      	end if

      	if(.not. (size(rho_pump,1) == N1_from_type(type_pump) .and. size(rho_pump,2) == N2_from_type(type_pump))) then
 			call print_error_message(-1,  'dimension error1 in integrate_pumped_rho')
 			stop
      	end if


		rho_pump = 0.0_dp
		rho_g(1,1) = 1
		call calculate_dipole_excitation('g','O',rho_g,rho_O,.true.)

		i0 = 1
      	do i = i0, Nt(1)
     		write(*,*) 'INTEGRUJI CAS T = ',(i-1)*dt*gt(1),' fs'

     		! since single_photon_I depends only on time difference, we can
     		! rho_pump(k,l,i) = rho_pump(k,l,i-1) + outer loop
     		if(i > i0) then
				do k=1,N1_from_type('E')
				do l=1,N2_from_type('E')
     				rho_pump(k,l,i) = rho_pump(k,l,i-1)
     			end do
     			end do
     		end if

      	do jj  = 0, i-i0
      	do jjj = i-jj-i0, 0, -1

      		if(.not.(jj == i-i0)) then
      			if(.not.(jjj == i-jj-i0)) then
      				exit
      			end if
      		end if

      		ttime  = (i-jj-1)*dt*gt(1)
      		tttime = (i-jj-jjj-1)*dt*gt(1)

      		if(single_photon_I_thermal(ttime,tttime,tau,w)/single_photon_I_thermal(ttime,ttime,tau,w) < 1e-2) then
     			cycle
    		end if

      		rho_O2 = 0.0_dp

			do k=1,N1_from_type('O')
			do l=1,N2_from_type('O')

			do r=1,N1_from_type('O')
			do s=1,N2_from_type('O')

				rho_O2(k,l) = rho_O2(k,l) + evops(1,1)%Ueg(k,l,r,s, jjj+1)*rho_O(r,s)

			end do
			end do

			end do
			end do

			call calculate_dipole_excitation('O','E',rho_O2,rho_E,.false.)

			do k=1,N1_from_type('E')
			do l=1,N2_from_type('E')

			do r=1,N1_from_type('E')
			do s=1,N2_from_type('E')

				rho_pump(k,l,i) = rho_pump(k,l,i) 	 					&
					+ evops(1,1)%Uee(k,l,r,s,jj+1)*rho_E(r,s)			&
					* single_photon_I_thermal(ttime,tttime,tau,w)	&
					* gt(1)*dt * gt(1)*dt

			end do
			end do

			end do
			end do

      	end do
      	end do
      	end do

      	do i = i0, Nt(1)
      		rho_pump(:,:,i) = rho_pump(:,:,i) + transpose(conjg(rho_pump(:,:,i)))
      	end do

	end subroutine integrate_pumped_rho_thermal



	subroutine debug()
		complex(dpc), dimension(N1_from_type('O'),N2_from_type('O')) :: rho_O
		complex(dpc), dimension(N1_from_type('E'),N2_from_type('E')) :: rho_E, eig1, eig2
		complex(dpc), dimension(N1_from_type('g'),N2_from_type('g')) :: rho_g

		complex(dpc), dimension(N1_from_type('E'),N2_from_type('E'), Nt(1)) :: rho_pump

		integer		:: k,l,a,b,i
      	character(len=10)	:: cha,chb
      	character(len=64)	:: name
      	character			:: type

      	real(dp) :: tau, w

		type = 'E'

		w 	=  0.0_dp
		tau	= 10.0_dp

		!toto uz snad funguje
		!call integrate_pumped_rho_thermal('E',rho_pump,tau,w)
		call integrate_pumped_rho_laser('E',rho_pump,tau,w,50.0_dp)

      	do a = 1, N1_from_type(type)
      	do b = 1, N2_from_type(type)

	      	write(cha,'(i1)') a
    	  	write(chb,'(i1)') b

			! v excitonove bazi
			if (type == 'E') then
				name = 'pump_E_'//trim(cha)//'_'//trim(chb)//'.dat'
      		else if (type == 'O') then
  	    		name = 'pump_O_'//trim(cha)//'_'//trim(chb)//'.dat'
    	  	else if (type == '2') then
      			name = 'pump_2_'//trim(cha)//'_'//trim(chb)//'.dat'
      		end if

			open(unit=11,file=trim(file_join(out_dir,trim(name))))

    	  	do i = 1, Nt(1)
				write(11,*) (i-1)*dt*gt(1), real(rho_pump(a,b,i)), aimag(rho_pump(a,b,i))
	      	end do

    	  	close(unit=11)


    	  	if(.not.(a == b)) then
    	  		cycle
    	  	end if


    	  	! v bazi vlastnich hodnot
			if (type == 'E') then
				name = 'pump_eig_E_'//trim(cha)//'_'//trim(chb)//'.dat'
      		else if (type == 'O') then
  	    		name = 'pump_eig_O_'//trim(cha)//'_'//trim(chb)//'.dat'
    	  	else if (type == '2') then
      			name = 'pump_eig_2_'//trim(cha)//'_'//trim(chb)//'.dat'
      		end if

!			open(unit=12,file=trim(file_join(out_dir,trim(name))))

!    	  	do i = 1, Nt(1)
!   	  		rho_E = rho_pump(:,:,i)
!				call spec(rho_E,eig1,eig2)
!				write(12,*) (i-1)*dt*gt(1), real(eig2(a,b)), aimag(eig2(a,b))
!	      	end do

!   	  	close(unit=12)

      	end do
      	end do


	end subroutine debug

end module module_qme

