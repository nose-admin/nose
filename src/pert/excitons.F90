#include "util_allocation.h"

! This modules prepares all variables associated with exciton model
! using inputs from wizards
!
module excitons

    use input_wizard

    implicit none

    real(dp), dimension(:), allocatable :: r

    ! energy shift due to static disorder
    real(dp), dimension(:,:), allocatable :: den

    private :: init_block_disorder, r, den

    character(len=64), private :: cbuff


    real(dp), dimension(:), allocatable :: ens

contains

    !
    ! Initializes excitons: makes sure that the current_e_block contains excitonic
    ! version of the current_s_block variables
    !
    ! Here we ONLY ALLOCATE. No calculations, except things that don't change between runs
    !
    subroutine init_excitons(err)
        integer, intent(out) :: err

        real(dp), dimension(:,:), allocatable :: HH, He
        real(dp), dimension(:,:), allocatable :: HH2, He2
        real(dp), dimension(:,:), pointer      :: SS => NULL(), S1 => NULL()
        real(dp), dimension(:,:), pointer      :: SS2 => NULL(), S12 => NULL()

        integer, dimension(:,:), pointer       :: itwo => NULL()
        integer, dimension(:), pointer         :: ione1 => NULL(), ione2 => NULL()

        real(dp), dimension(:), pointer        :: e => NULL(), e2 => NULL(), eg => NULL()
        real(dp), dimension(:,:), pointer      :: dx => NULL(),dy => NULL(),dz => NULL(),dd => NULL()
        real(dp), dimension(:,:), pointer      :: dx2 => NULL(),dy2 => NULL(),dz2 => NULL(),dd2 => NULL()
        real(dp), dimension(3)                  :: vp
        real(dp), dimension(:), allocatable    :: flf
        real(dp), dimension(:), pointer        :: fl_fac => NULL()
        real(dp) :: sum


		real(dp), dimension(:,:,:), allocatable :: k_bands

        integer :: i,j,m,n,k,l, bk, nU, mU

        integer :: N2, Nkk

        integer :: N_tot, bcS, bcE

        err = 0

		call print_log_message("Exciton module ...", 8)

        bcS = get_byte_count()

        !
        ! initialization of the disorder
        !
        call init_random(rseed)

        !
        ! Allocated total fluorescence factors
        !
        N_tot = Block_Sizes_Added(nr_blocks)
!        if (.not.allocated(flf)) then
!            !allocate(flf(N_tot))
!            ALLOCATE(flf,(N_tot))
!            ! register_allocation(flf)
!        end if

        !
        ! set blocks to the start, use site representation
        !
        call resources_rewind_blocks()
        call resources_set_all_pointers(RESOURCES_SITE_REP)


		if ( current_s_block%periodic == 1 ) then

			Nkk = current_s_block%periodicity


			if (.not.allocated(HH)) then
				ALLOCATE(HH,(N1,N1))
				ALLOCATE(SS,(N1,N1))
				ALLOCATE(S1,(N1,N1))
				ALLOCATE(He,(N1,N1))
			end if

			ALLOCATE(Bands_1ex,(N1,Nkk,Nkk))
			do m = 1, Nkk
				do n = 1, Nkk
					nU = n - Nkk/2
					mU = m - Nkk/2
					HH = 0.0_dp
					! one exciton part
    		        do i = 1, N1
            		    HH(i,i) = current_s_block%en_orig(i) ! + energy_disorder(i)
                		do j = 1, N1
                			HH(i,j) = HH(i,j) + current_s_block%Jcell(1,i,j)*2*cos(TWOPI_D*mU/Nkk)
                			HH(i,j) = HH(i,j) + current_s_block%Jcell(2,i,j)*2*cos(TWOPI_D*nU/Nkk)
                    		HH(i,j) = HH(i,j) + current_s_block%Jcell(3,i,j)*2*cos(TWOPI_D*(nU+mU)/Nkk)
                    		HH(i,j) = HH(i,j) + current_s_block%Jcell(4,i,j)*2*cos(TWOPI_D*(nU-mU)/Nkk)
                    		if (i/=j) then
                        		HH(i,j) = HH(i,j) + current_s_block%J(i,j) !+ coupling_disorder(i)
							end if
                		end do
            		end do

            		call spec(HH,SS,He)
            		call inv(SS,S1)

					do k = 1, N1
						Bands_1ex(k,m,n) = He(k,k)
					end do

				end do
			end do

   		    !
        	! some things will be logged (only master process)
        	!
        	if ( (parallel_id == 0).and.(i_realization == 1) ) then
            	open(unit=11,file= trim(file_join("log","bands_1ex_1.dat")),status="replace")
            	write(11,'(a)') "# One exciton bands: "
            	write(11,'(a)') "# "
				do i = 1, Nkk
				!do j = 1, Nkk
					write(11,*) Bands_1ex(1,i,:)
				!end do
				end do
            	open(unit=12,file=trim(file_join("log","bands_1ex_2.dat")))
            	write(12,'(a)') "# Two-exciton bands: "
            	write(12,'(a)') "# "
				do i = 1, Nkk
				!do j = 1, Nkk
					write(12,*) Bands_1ex(2,i,:)
				!end do
				end do
            	close(11)
				close(12)

        	end if




			stop

		else
        !
        !
        !
        do

            if (current_s_block%id == 0) exit

            !
            ! create new exciton block for each site block
            !
            call resources_new_ex_block()

            if (.not.resources_have_next_block()) exit

            call resources_next_block()
            call resources_set_all_pointers(RESOURCES_SITE_REP)

        end do


        !
        ! reset everything to the beginning
        !
        call resources_rewind_blocks()
        call resources_set_all_pointers(RESOURCES_SITE_REP)


        !
        ! loop over blocks and set excitonic properties
        !
        bk = 0
        do

            bk = bk + 1

            if (current_s_block%id == 0) exit

            !
            ! initialize disorder in current block
            !
            call init_block_disorder(current_s_block)

            ! number of 2 exciton states
            N2 = N1*(N1-1)/2

            if (.not.allocated(HH)) then
                if (use_twoexcitons) then
                    ALLOCATE_S(current_e_block%N2)
                    current_e_block%N2 = N2
                end if
            end if
            ALLOCATE(SS,(N1,N1))
            ALLOCATE(S1,(N1,N1))
            ALLOCATE(e,(N1))
            ALLOCATE(eg,(1))
            if (use_twoexcitons) then
                ALLOCATE(SS2,(N2,N2))
                ALLOCATE(S12,(N2,N2))
                ALLOCATE(e2,(N2))
            end if

            ALLOCATE(current_e_block%ll,(N1))


            ALLOCATE(current_e_block%gg,(N1,Nte))
            ALLOCATE(current_e_block%cc,(N1,Nte))

            ! two exciton part
            if (use_twoexcitons) then

                !
                ! initialize conversion matrices
                !
                !call allocate_pointer(itwo,N1,N1)
                !call allocate_pointer(ione1,N1)
                !call allocate_pointer(ione2,N1)
                ALLOCATE(itwo,(N1,N1))
                ALLOCATE(ione1,(N1*(N1-1)/2))
                ALLOCATE(ione2,(N1*(N1-1)/2))
                k = 0
                itwo = 0
                ione1 = 0
                ione2 = 0
                do i = 1, N1
                    do j = i+1, N1
                        k = k + 1
                        itwo(i,j) = k
                        ione1(k) = i
                        ione2(k) = j
                    end do
                end do
                current_e_block%itwo => itwo
                current_e_block%ione1 => ione1
                current_e_block%ione2 => ione2

                !
                ! create Hamiltonian
                !
 !               do i = 1, N2
 !                   k = ione1(i)
 !                   l = ione2(i)
 !                   HH2(i,i) = HH(k,k) + HH(l,l)
 !                   do j = 1, N2
 !                      if (i /= j) then
 !                        m = ione1(j)
 !                        n = ione2(j)
 !                        if (k==m) then
 !                           HH2(i,j) = HH(l,n)
 !                        else if (l==n) then
 !                           HH2(i,j) = HH(k,m)
 !                        else if (k==n) then
 !                           HH2(i,j) = HH(l,m)
 !                        else if (l==m) then
 !                           HH2(i,j) = HH(k,n)
 !                        end if
 !                      end if
 !                   end do
 !               end do

            end if

            !
            ! find eigenvalues and eigenvectors (transformation matrix)
            !
            !print *, N1, size(HH,1), size(HH,2), size(SS,1), size(SS,2), size(He,1), size(He,2)
 !           call spec(HH,SS,He)
 !           call inv(SS,S1)
 !           if (use_twoexcitons) then
 !               call spec(HH2,SS2,He2)
 !               call inv(SS2,S12)
 !           end if

            !
            ! save eigenvalues and transformation
            !
            current_e_block%SS => SS
            current_e_block%S1 => S1
!            do i = 1, N1
!                 e(i) = He(i,i)
!            end do
            if (use_twoexcitons) then
                current_e_block%SS_2 => SS2
                current_e_block%S1_2 => S12
!                do i = 1, N2
!                    e2(i) = He2(i,i)
!                end do
            end if



            !
            ! save eigenenergies
            !
            current_e_block%en => e
            current_e_block%eng => eg
!            do i = 1, N1
!                current_s_block%en(i) = HH(i,i)
!            end do
!            DEALLOCATE(HH)
!            DEALLOCATE(He)
            !call deregister_allocation(HH)
            !call deregister_allocation(He)
            !deallocate(HH,He)

            if (use_twoexcitons) then
                current_e_block%en_2 => e2
                !call deregister_allocation(HH2)
                !call deregister_allocation(He2)
                !deallocate(HH2,He2)
!                DEALLOCATE(HH2)
!                DEALLOCATE(He2)
            end if


            !
            ! calculate and save dipole moments
            !
            !call allocate_pointer(dx,N1)
            ALLOCATE(dx,(N1,1))
 !           dx = matmul(S1,current_s_block%dx)
            current_e_block%dx => dx
            !call allocate_pointer(dy,N1)
            ALLOCATE(dy,(N1,1))
 !           dy = matmul(S1,current_s_block%dy)
            current_e_block%dy => dy
            !call allocate_pointer(dz,N1)
            ALLOCATE(dz,(N1,1))
 !           dz = matmul(S1,current_s_block%dz)
            current_e_block%dz => dz

            if (use_twoexcitons) then
                !call allocate_pointer(dx2,N2,N1)
                !call allocate_pointer(dy2,N2,N1)
                !call allocate_pointer(dz2,N2,N1)
                ALLOCATE(dx2,(N2,N1))
                ALLOCATE(dy2,(N2,N1))
                ALLOCATE(dz2,(N2,N1))

!                do i = 1, N2
!                    k = ione1(i)
!                    l = ione2(i)
!                    do j = 1, N1
!                        if (j == k) then
!                            dx2(i,j) = current_s_block%dx(l)
!                            dy2(i,j) = current_s_block%dy(l)
!                            dz2(i,j) = current_s_block%dz(l)
!                        else if (j == l) then
!                            dx2(i,j) = current_s_block%dx(k)
!                            dy2(i,j) = current_s_block%dy(k)
!                            dz2(i,j) = current_s_block%dz(k)
!                        else
!                            dx2(i,j) = 0.0_dp
!                            dy2(i,j) = 0.0_dp
!                            dz2(i,j) = 0.0_dp
!                        end if
!                    end do
!                end do

                ! transform it into excitonic representation

 !               dx2 = matmul(S12,matmul(dx2,SS))
 !               dy2 = matmul(S12,matmul(dy2,SS))
 !               dz2 = matmul(S12,matmul(dz2,SS))

                current_e_block%dx_2 => dx2
                current_e_block%dy_2 => dy2
                current_e_block%dz_2 => dz2

            end if

            !
            ! calculate save dipole moment lengths
            !
            !call allocate_pointer(dd,N1)
            ALLOCATE(dd,(N1,1))
            ! find spatial average
!            do i = 1, N1
!                dd(i) = sqrt(current_e_block%dx(i)**2 + &
!                current_e_block%dy(i)**2 + current_e_block%dz(i)**2)
!               ! print *, current_e_block%en(i), dd(i)
!            end do

            current_e_block%dd => dd

            if (use_twoexcitons) then
                !call allocate_pointer(dd2,N2,N1)
                ALLOCATE(dd2,(N2,N1))
!                do i = 1, N2
!                    do j = 1, N1
!                        dd2(i,j) = sqrt(current_e_block%dx_2(i,j)**2 + &
!                        current_e_block%dy_2(i,j)**2 + current_e_block%dz_2(i,j)**2)
!                        !print *, i,j, dd2(i,j)
!                    end do
!                   ! print *, current_e_block%en(i), dd(i)
!                end do
                current_e_block%dd_2 => dd2
            end if

!            if ((parallel_id == 0).and.(i_realization == 1)) then
!
!                write(11,'(a,i5)') "# Block nr. ", current_e_block%id
!                write(11,'(a)') "#"
!                write(11,'(a)') "# Transition energy [cm^-1]      Dipole strength [D] "
!                write(11,'(a)') "# -----------------------------------------------------"
!                do i = 1,N1
!
!                       write(11,'(f18.6,e18.6)') &
!                          (current_e_block%en(i))*Energy_internal_to_cm, &
!                          current_e_block%dd(i)**2 !1.76*(current_e_block%dd(i)**2)/(current_e_block%dd(1)**2)
!
!                end do
!                write(11,'(a)') "# "
!                write(11,'(a)') "# "
!
!                if (use_twoexcitons) then
!                write(12,'(a,i5)') "# Block nr. ", current_e_block%id
!                write(12,'(a)') "#"
!                write(12,'(a)') "# Transition energy [cm^-1]      Dipole strength [D] "
!                write(12,'(a)') "# -----------------------------------------------------"
!                do i = 1,N2
!                    do j = 1, N1
!                       write(12,'(2i5,e18.6,e18.6)') i, j, &
!                          (current_e_block%en_2(i)-current_e_block%en(j))*Energy_internal_to_cm, &
!                          current_e_block%dd_2(i,j)**2 !1.76*(current_e_block%dd(i)**2)/(current_e_block%dd(1)**2)
!                    end do
!                end do
!                write(12,'(a)') "# "
!                write(12,'(a)') "# "
!                end if
!
!
!            end if


            !
            ! rotational strength
            !


            !call allocate_pointer(rm,N1)
            ALLOCATE(rm,(N1))
!            rm = matmul(S1,current_s_block%rs)

            ! CONSERVATIVE PART OF THE ROTATIONAL STREGTH
!            do i = 1, N1
!                do n = 1, N1
!                    do m = n+1,N1
!                        vp(1) = dy(m)*dz(n) - dz(m)*dy(n)
!                        vp(2) = dz(m)*dx(n) - dx(m)*dz(n)
!                        vp(3) = dx(m)*dy(n) - dy(m)*dx(n)
!                        rm(i) = rm(i) + S1(i,m)*S1(i,n)*dot_product( &
!                         current_s_block%rr(m,:)-current_s_block%rr(n,:),vp)
!                    end do
!                end do
!            end do


            !
            ! save rotational strength
            !
            current_e_block%rm => rm

!            flf(Block_Sizes_Previous(bk)+1:Block_Sizes_Added(bk)) = current_e_block%en(1:N1)

            if (.not.resources_have_next_block()) exit

            call resources_next_block()
            call resources_set_all_pointers(RESOURCES_SITE_REP)

        end do

        !
        ! Close logging
        !
!        if ((parallel_id == 0).and.(i_realization == 1)) then
!
!            close(unit=11)
!            close(unit=12)
!
!        end if

        call random_get_seed(rseed)


        !
        ! Create fuorescence population factors
        !
!        do i = 1, N_tot
!           sum = sum + exp(-flf(i)/(kB_cmK*Energy_cm_to_internal*temp))
!        end do

        call resources_rewind_blocks()
        call resources_set_all_pointers(RESOURCES_SITE_REP)

        k = 0
        do

            if (current_s_block%id == 0) exit

            !call allocate_pointer(fl_fac,N1)
            ALLOCATE(fl_fac,(N1))

!            do i = 1, N1
!                k = k + 1
!                fl_fac(i) = exp(-flf(k)/(kB_cmK*Energy_cm_to_internal*temp))/sum
!            end do

            !
            ! save fluorescence population factors
            !
            current_e_block%fl_fac => fl_fac

            if (.not.resources_have_next_block()) exit

            call resources_next_block()
            call resources_set_all_pointers(RESOURCES_SITE_REP)

       end do

       bcE = get_byte_count()

!       write(cbuff,'(a,i8,a)') "Total of ",bcE-bcS," bytes allocated while initializing excitons"
!       call print_log_message(trim(cbuff),5)
!       write(cbuff,'(a,f8.4,a)') "Currently using ", real(get_byte_count())/(1024.0_dp**2), " MB"
!       call print_log_message(trim(cbuff),5)

	end if

    call print_log_message("... exciton module done", 8)


	end subroutine init_excitons



    !
    ! Initializes excitons: makes sure that the current_e_block contains excitonic
    ! version of the current_s_block variables
    !
    subroutine update_excitons(err)
        integer, intent(out) :: err

        real(dp), dimension(:,:), allocatable :: HH, He
        real(dp), dimension(:,:), allocatable :: HH2, He2
        real(dp), dimension(:,:), allocatable      :: SS, S1
        real(dp), dimension(:,:), allocatable      :: SS2, S12

        integer, dimension(:,:), pointer       :: itwo => NULL()
        integer, dimension(:), pointer         :: ione1 => NULL(), ione2 => NULL()

        real(dp), dimension(:), pointer        :: e => NULL(), e2 => NULL()
        real(dp), dimension(:,:), pointer      :: dx => NULL(),dy => NULL(),dz => NULL(),dd => NULL()
        real(dp), dimension(:,:), pointer      :: dx2 => NULL(),dy2 => NULL(),dz2 => NULL(),dd2 => NULL()
        real(dp), dimension(3)                  :: vp
        real(dp), dimension(:), allocatable    :: flf
        real(dp), dimension(:), pointer        :: fl_fac => NULL()
        real(dp) :: sum, jcp


		real(dp), dimension(:,:,:), allocatable :: k_bands

        integer :: i,j,m,n,k,l, bk, nU, mU

        integer :: N2, Nkk

        integer :: N_tot, bcS, bcE

        err = 0

        bcS = get_byte_count()

		call print_log_message("Update excitons ...", 8)

        !
        ! initialization of the disorder
        !
        call init_random(rseed)

        !
        ! Allocated total fluorescence factors
        !
        N_tot = Block_Sizes_Added(nr_blocks)
        if (.not.allocated(flf)) then
            !allocate(flf(N_tot))
            ALLOCATE(flf,(N_tot))
            ! register_allocation(flf)
        end if

        !
        ! set blocks to the start, use site representation
        !
        call resources_rewind_blocks()
        call resources_set_all_pointers(RESOURCES_SITE_REP)


		if ( current_s_block%periodic == 1 ) then

			Nkk = current_s_block%periodicity


			if (.not.allocated(HH)) then
				ALLOCATE(HH,(N1,N1))
				ALLOCATE(SS,(N1,N1))
				ALLOCATE(S1,(N1,N1))
				ALLOCATE(He,(N1,N1))
			end if

			ALLOCATE(Bands_1ex,(N1,Nkk,Nkk))
			do m = 1, Nkk
				do n = 1, Nkk
					nU = n - Nkk/2
					mU = m - Nkk/2
					HH = 0.0_dp
					! one exciton part
    		        do i = 1, N1
            		    HH(i,i) = current_s_block%en_orig(i) ! + energy_disorder(i)
                		do j = 1, N1
                			HH(i,j) = HH(i,j) + current_s_block%Jcell(1,i,j)*2*cos(TWOPI_D*mU/Nkk)
                			HH(i,j) = HH(i,j) + current_s_block%Jcell(2,i,j)*2*cos(TWOPI_D*nU/Nkk)
                    		HH(i,j) = HH(i,j) + current_s_block%Jcell(3,i,j)*2*cos(TWOPI_D*(nU+mU)/Nkk)
                    		HH(i,j) = HH(i,j) + current_s_block%Jcell(4,i,j)*2*cos(TWOPI_D*(nU-mU)/Nkk)
                    		if (i/=j) then
                        		HH(i,j) = HH(i,j) + current_s_block%J(i,j) !+ coupling_disorder(i)
							end if
                		end do
            		end do

            		call spec(HH,SS,He)
            		call inv(SS,S1)

					do k = 1, N1
						Bands_1ex(k,m,n) = He(k,k)
					end do

				end do
			end do

   		    !
        	! some things will be logged (only master process)
        	!
        	if ( (parallel_id == 0).and.(i_realization == 1) ) then
            	open(unit=11,file= trim(file_join("log","bands_1ex_1.dat")),status="replace")
            	write(11,'(a)') "# One exciton bands: "
            	write(11,'(a)') "# "
				do i = 1, Nkk
				!do j = 1, Nkk
					write(11,*) Bands_1ex(1,i,:)
				!end do
				end do
            	open(unit=12,file=trim(file_join("log","bands_1ex_2.dat")))
            	write(12,'(a)') "# Two-exciton bands: "
            	write(12,'(a)') "# "
				do i = 1, Nkk
				!do j = 1, Nkk
					write(12,*) Bands_1ex(2,i,:)
				!end do
				end do
            	close(11)
				close(12)

        	end if




			stop

		else
        !
        !
        !


        !
        ! reset everything to the beginning
        !
        call resources_rewind_blocks()
        call resources_set_all_pointers(RESOURCES_SITE_REP)

        !
        ! some things will be logged (only master process)
        !
        if ((parallel_id == 0).and.(i_realization == 1)) then
            open(unit=11,file= trim(file_join("log","transitions_01.dat")),status="replace")
            write(11,'(a)') "# Optical transitions: "
            write(11,'(a)') "# "
            open(unit=12,file=trim(file_join("log","transitions_12.dat")))
            write(12,'(a)') "# Optical transitions: "
            write(12,'(a)') "# "
        end if


        !
        ! loop over blocks and set excitonic properties
        !
        bk = 0
        do

            bk = bk + 1

            if (current_s_block%id == 0) exit


            !
            ! initialize disorder in current block
            !
            call init_block_disorder(current_s_block)

            ! number of 2 exciton states
            N2 = N1*(N1-1)/2

            if (.not.allocated(HH)) then
                !allocate(HH(N1,N1),He(N1,N1))
                ALLOCATE(HH,(N1,N1))
                ALLOCATE(He,(N1,N1))
                !call register_allocation(HH)
                !call register_allocation(He)
                if (use_twoexcitons) then
                    !allocate(HH2(N2,N2),He2(N2,N2))
                    ALLOCATE(HH2,(N2,N2))
                    ALLOCATE(He2,(N2,N2))
!! done in init
!                    ALLOCATE_S(current_e_block%N2)
!                    current_e_block%N2 = N2
                    !call register_allocation(HH2)
                    !call register_allocation(He2)
                end if
            end if
            !call allocate_pointer(SS,N1,N1)
            !call allocate_pointer(S1,N1,N1)
            !call allocate_pointer(e,N1)
            ALLOCATE(SS,(N1,N1))
            ALLOCATE(S1,(N1,N1))
            ALLOCATE(e,(N1))
            if (use_twoexcitons) then
                !call allocate_pointer(SS2,N2,N2)
                !call allocate_pointer(S12,N2,N2)
                !call allocate_pointer(e2,N1)
                ALLOCATE(SS2,(N2,N2))
                ALLOCATE(S12,(N2,N2))
                ALLOCATE(e2,(N2))
            end if

            HH = 0.0_dp
            SS = 0.0_dp
            S1 = 0.0_dp
            if (use_twoexcitons) then
                HH2 = 0.0_dp
                SS2 = 0.0_dp
                S12 = 0.0_dp
            end if


            !
            ! make block Hamiltonian
            !

            ! one exciton part
            do i = 1, N1
                HH(i,i) = current_s_block%en_orig(i) + energy_disorder(i)
                do j = 1, N1
                    if (i/=j) then
                    	jcp = coupling_disorder(i,j)
                    	if (jj(i,j) < 0.0d0) then
                        	HH(i,j) = jj(i,j) - jcp
                        else
                        	HH(i,j) = jj(i,j) + jcp
                        end if
                    end if
                end do
            end do

            ! two exciton part
            if (use_twoexcitons) then

                ALLOCATE(itwo,(N1,N1))
                ALLOCATE(ione1,(N1*(N1-1)/2))
                ALLOCATE(ione2,(N1*(N1-1)/2))
!                k = 0
                itwo = current_e_block%itwo
                ione1 = current_e_block%ione1
                ione2 = current_e_block%ione2

                !
                ! create Hamiltonian
                !
                do i = 1, N2
                    k = ione1(i)
                    l = ione2(i)
                    HH2(i,i) = HH(k,k) + HH(l,l)
                    do j = 1, N2
                       if (i /= j) then
                         m = ione1(j)
                         n = ione2(j)
                         if (k==m) then
                            HH2(i,j) = HH(l,n)
                         else if (l==n) then
                            HH2(i,j) = HH(k,m)
                         else if (k==n) then
                            HH2(i,j) = HH(l,m)
                         else if (l==m) then
                            HH2(i,j) = HH(k,n)
                         end if
                       end if
                    end do
                end do

            end if

            !
            ! find eigenvalues and eigenvectors (transformation matrix)
            !
            !print *, N1, size(HH,1), size(HH,2), size(SS,1), size(SS,2), size(He,1), size(He,2)
            call spec(HH,SS,He)

			!print *, He*Energy_internal_to_cm;

            call inv(SS,S1)
            if (use_twoexcitons) then
                call spec(HH2,SS2,He2)
                call inv(SS2,S12)
            end if

            !
            ! save eigenvalues and transformation
            !
            current_e_block%SS = SS
            current_e_block%S1 = S1
            do i = 1, N1
                 e(i) = He(i,i)
            end do
            if (use_twoexcitons) then
                current_e_block%SS_2 = SS2
                current_e_block%S1_2 = S12
                do i = 1, N2
                    e2(i) = He2(i,i)
                end do
            end if


            !
            ! save eigenenergies
            !
            current_e_block%en = e
            current_e_block%eng = 0.0_dp
            do i = 1, N1
                current_s_block%en(i) = HH(i,i)
            end do
            DEALLOCATE(HH)
            DEALLOCATE(He)
            DEALLOCATE(e)
            !call deregister_allocation(HH)
            !call deregister_allocation(He)
            !deallocate(HH,He)

            if (use_twoexcitons) then
                current_e_block%en_2 = e2
                !call deregister_allocation(HH2)
                !call deregister_allocation(He2)
                !deallocate(HH2,He2)
                DEALLOCATE(HH2)
                DEALLOCATE(He2)
                DEALLOCATE(e2)
            end if


            !
            ! calculate and save dipole moments
            !
            !call allocate_pointer(dx,N1)
            ALLOCATE(dx,(N1,1))
            dx = matmul(S1,current_s_block%dx)
            current_e_block%dx = dx
            !call allocate_pointer(dy,N1)
            ALLOCATE(dy,(N1,1))
            dy = matmul(S1,current_s_block%dy)
            current_e_block%dy = dy
            !call allocate_pointer(dz,N1)
            ALLOCATE(dz,(N1,1))
            dz = matmul(S1,current_s_block%dz)
            current_e_block%dz = dz

            if (use_twoexcitons) then

                ALLOCATE(dx2,(N2,N1))
                ALLOCATE(dy2,(N2,N1))
                ALLOCATE(dz2,(N2,N1))

                do i = 1, N2
                    k = ione1(i)
                    l = ione2(i)
                    do j = 1, N1
                        if (j == k) then
                            dx2(i,j) = current_s_block%dx(l,1)
                            dy2(i,j) = current_s_block%dy(l,1)
                            dz2(i,j) = current_s_block%dz(l,1)
                        else if (j == l) then
                            dx2(i,j) = current_s_block%dx(k,1)
                            dy2(i,j) = current_s_block%dy(k,1)
                            dz2(i,j) = current_s_block%dz(k,1)
                        else
                            dx2(i,j) = 0.0_dp
                            dy2(i,j) = 0.0_dp
                            dz2(i,j) = 0.0_dp
                        end if
                    end do
                end do

                ! transform it into excitonic representation

                dx2 = matmul(S12,matmul(dx2,SS))
                dy2 = matmul(S12,matmul(dy2,SS))
                dz2 = matmul(S12,matmul(dz2,SS))

                current_e_block%dx_2 = dx2
                current_e_block%dy_2 = dy2
                current_e_block%dz_2 = dz2

            end if


            !
            ! calculate save dipole moment lengths
            !
            !call allocate_pointer(dd,N1)
            ALLOCATE(dd,(N1,1))
            ! find spatial average
            do i = 1, N1
                dd(i,1) = sqrt(current_e_block%dx(i,1)**2 + &
                current_e_block%dy(i,1)**2 + current_e_block%dz(i,1)**2)
               ! print *, current_e_block%en(i), dd(i)
            end do

            current_e_block%dd = dd

            if (use_twoexcitons) then
                !call allocate_pointer(dd2,N2,N1)
                ALLOCATE(dd2,(N2,N1))
                do i = 1, N2
                    do j = 1, N1
                        dd2(i,j) = sqrt(current_e_block%dx_2(i,j)**2 + &
                        current_e_block%dy_2(i,j)**2 + current_e_block%dz_2(i,j)**2)
                        !print *, i,j, dd2(i,j)
                    end do
                   ! print *, current_e_block%en(i), dd(i)
                end do
                current_e_block%dd_2 = dd2
            end if

            if ((parallel_id == 0).and.(i_realization == 1)) then

                write(11,'(a,i5)') "# Block nr. ", current_e_block%id
                write(11,'(a)') "#"
                write(11,'(a)') "# Transition energy [cm^-1]      Dipole strength [D] "
                write(11,'(a)') "# -----------------------------------------------------"
                do i = 1,N1

                       write(11,'(f18.6,e18.6)') &
                          (current_e_block%en(i))*Energy_internal_to_cm, &
                          current_e_block%dd(i,1)**2 !1.76*(current_e_block%dd(i)**2)/(current_e_block%dd(1)**2)

                end do
                write(11,'(a)') "# "
                write(11,'(a)') "# "

                if (use_twoexcitons) then
                write(12,'(a,i5)') "# Block nr. ", current_e_block%id
                write(12,'(a)') "#"
                write(12,'(a)') "# Transition energy [cm^-1]      Dipole strength [D] "
                write(12,'(a)') "# -----------------------------------------------------"
                do i = 1,N2
                    do j = 1, N1
                       write(12,'(2i5,e18.6,e18.6)') i, j, &
                          (current_e_block%en_2(i)-current_e_block%en(j))*Energy_internal_to_cm, &
                          current_e_block%dd_2(i,j)**2 !1.76*(current_e_block%dd(i)**2)/(current_e_block%dd(1)**2)
                    end do
                end do
                write(12,'(a)') "# "
                write(12,'(a)') "# "
                end if


            end if


            !
            ! rotational strength
            !


            !call allocate_pointer(rm,N1)
            ALLOCATE(rm,(N1))
            rm = matmul(S1,current_s_block%rs)

            ! CONSERVATIVE PART OF THE ROTATIONAL STREGTH
            do i = 1, N1
                do n = 1, N1
                    do m = n+1,N1
                        vp(1) = dy(m,1)*dz(n,1) - dz(m,1)*dy(n,1)
                        vp(2) = dz(m,1)*dx(n,1) - dx(m,1)*dz(n,1)
                        vp(3) = dx(m,1)*dy(n,1) - dy(m,1)*dx(n,1)
                        rm(i) = rm(i) + S1(i,m)*S1(i,n)*dot_product( &
                         current_s_block%rr(m,:)-current_s_block%rr(n,:),vp)
                    end do
                end do
            end do

            !
            ! save rotational strength
            !
            current_e_block%rm = rm

            DEALLOCATE(SS)
            DEALLOCATE(S1)
			DEALLOCATE(dx)
			DEALLOCATE(dy)
			DEALLOCATE(dz)
			DEALLOCATE(dd)
			DEALLOCATE(rm)
			if (use_twoexcitons) then
                DEALLOCATE(SS2)
                DEALLOCATE(S12)
				DEALLOCATE(itwo)
				DEALLOCATE(ione1)
				DEALLOCATE(ione2)
				DEALLOCATE(dx2)
				DEALLOCATE(dy2)
				DEALLOCATE(dz2)
				DEALLOCATE(dd2)
			end if

            flf(Block_Sizes_Previous(bk)+1:Block_Sizes_Added(bk)) = current_e_block%en(1:N1)

            if (.not.resources_have_next_block()) exit

            call resources_next_block()
            call resources_set_all_pointers(RESOURCES_SITE_REP)

        end do

        !
        ! Close logging
        !
        if ((parallel_id == 0).and.(i_realization == 1)) then

            close(unit=11)
            close(unit=12)

        end if

        call random_get_seed(rseed)


        !
        ! Create fuorescence population factors
        !
        do i = 1, N_tot
           sum = sum + exp(-flf(i)/(kB_cmK*Energy_cm_to_internal*temp))
        end do

        call resources_rewind_blocks()
        call resources_set_all_pointers(RESOURCES_SITE_REP)

        k = 0
        do

            if (current_s_block%id == 0) exit

            !call allocate_pointer(fl_fac,N1)
            ALLOCATE(fl_fac,(N1))

            do i = 1, N1
                k = k + 1
                fl_fac(i) = exp(-flf(k)/(kB_cmK*Energy_cm_to_internal*temp))/sum
            end do

            !
            ! save fluorescence population factors
            !
            current_e_block%fl_fac = fl_fac

            if (.not.resources_have_next_block()) exit

            call resources_next_block()
            call resources_set_all_pointers(RESOURCES_SITE_REP)

       end do

       DEALLOCATE(flf)
	   DEALLOCATE(fl_fac)

       bcE = get_byte_count()

!       write(cbuff,'(a,i8,a)') "Total of ",bcE-bcS," bytes allocated while initializing excitons"
!       call print_log_message(trim(cbuff),5)
!       write(cbuff,'(a,f8.4,a)') "Currently using ", real(get_byte_count())/(1024.0_dp**2), " MB"
!       call print_log_message(trim(cbuff),5)

	end if


	call print_log_message("... exciton upate done", 8)




    end subroutine update_excitons




    !
    ! Creates energy shifts due to the disorder
    !
    subroutine init_block_disorder(csb)
        type(site_block), pointer :: csb
        integer :: i,j

        if (allocated(r)) then
             DEALLOCATE(r)
        end if
        ALLOCATE(r,(csb%N1**2))

        call random_Normal(r)

        if (allocated(den)) then
            DEALLOCATE(den)
        end if
        ALLOCATE(den,(csb%N1,csb%N1))

		! first realization is always without any disorder
        if ((parallel_id == 0).and.(i_realization == 1)) then

            den = 0.0_dp

        else

        do i = 1, csb%N1
        do j = i, csb%N1
        	if (i == j) then
        		den(i,i) = csb%dwidth(i,i)*r( (i-1)*csb%N1 + i )
        	else

        		! here we fix one r for each realization --- corresponds to perfectly correlated disorder
            	den(i,j) = csb%dwidth(i,j)*r( (1-1)*csb%N1 + 2 )
            	den(j,i) = den(i,j)
           	end if
        end do
        end do



        end if

    end subroutine init_block_disorder

    !
    !
    !
    function energy_disorder(i) result (res)
        integer, intent(in) :: i
        real(dp) :: res
        res = den(i,i)
    end function energy_disorder

    !
    !
    !
    function coupling_disorder(i,j) result (res)
        integer, intent(in) :: i,j
        real(dp) :: res
        res = den(i,j) !0.0_dp
    end function coupling_disorder


end module excitons
