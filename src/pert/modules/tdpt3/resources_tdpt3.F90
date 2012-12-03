#include "util_allocation.h"
!
!
! Common resources for the TDPT-3 module
!
!
module resources_tdpt3

    use resources
    implicit none

    complex(dpc), dimension(:,:), pointer     :: Ueg ! coherence evolution operator U_{eg}(t)


    type e_ops_storage
        integer                               :: id
        complex(dpc), dimension(:,:), pointer :: Ueg
        complex(dpc), dimension(:,:), pointer :: feg
        type(e_ops_storage), pointer          :: next
    end type

    type(e_ops_storage), pointer              :: storage
    type(e_ops_storage), pointer              :: current_Ueg

	complex(dpc), dimension(:,:,:), pointer :: Ufe      ! coherence evolution operator U_{fe}(t)

 !PGF95   private                                   :: e_ops_storage
 !PGF95   private                                   :: storage

    logical                                   :: polar_1_collected
    logical                                   :: polar_fl_collected
    logical                                   :: tdep_cd_collected
    logical                                   :: hist_relax_collected
    logical                                   :: pump_probe_collected

contains

    !
    ! Initialization of the module resources
    !
    subroutine init_resources_tdpt3()


        call allocate_evolution_operators()

        polar_1_collected = .false.
        tdep_cd_collected = .false.
        polar_fl_collected = .false.
        hist_relax_collected = .false.
        pump_probe_collected = .false.


    end subroutine init_resources_tdpt3

    !
    ! Initialization of the module resources
    !
    subroutine update_resources_tdpt3()


        call prepare_evolution_operators()

        polar_1_collected = .false.
        tdep_cd_collected = .false.
        polar_fl_collected = .false.
        hist_relax_collected = .false.
        pump_probe_collected = .false.


    end subroutine update_resources_tdpt3


    subroutine allocate_evolution_operators

        complex(dpc), dimension(:,:), pointer :: lUeg, fUeg
        integer :: i,k,t, kk
        real(dp) :: tt
        real(dp), dimension(:), pointer :: dr
        integer :: N2


        ! ***********************************************************

        allocate(storage)
       	storage%id = 0
        storage%next => null()
        storage%Ueg  => null()
        storage%feg  => null()


        call resources_rewind_blocks()
        call resources_set_all_pointers(RESOURCES_EXCITON_REP)
        call rewind_evolution_operators()


        call new_evolution_operators()

        do

            ALLOCATE(lUeg,(N1,Nt(1)))
            ALLOCATE(fUeg,(N1,Nt(1)))

            current_Ueg%Ueg => lUeg
            current_Ueg%feg => fUeg

            if (.not.resources_have_next_block()) exit

            call resources_next_block()
            call resources_set_all_pointers(RESOURCES_EXCITON_REP)
!            call new_evolution_operators()


        end do

        ! ***************************************************************************


        ! new scheme of storing evolution operator
        call resources_rewind_blocks()
        call resources_set_all_pointers(RESOURCES_EXCITON_REP)
        !allocate(dr(N1))

        i = 1

        do



            N2 = N1*(N1-1)/2
            ALLOCATE(evops(i,i)%Ueg,(N1, 1,N1, 1,Nt(1)))
            ALLOCATE(evops(i,i)%Uee,(N1,N1,N1,N1,Nt(2)))
!*            ALLOCATE(evops(i,i)%Ufe,(N2,N1,N2,N1,Nt(1)))
            ALLOCATE(evops(i,i)%UfeS,(N2,N1,Nt(1)))
            ALLOCATE(evops(i,i)%Ugg,(1,1,1,1,Nt(1)))
            evops(i,i)%Ueg = 0.0_dp
            call super_unity(N1,Nt(2),evops(i,i)%Uee)
!*            evops(i,i)%Ufe = 0.0_dp
			evops(i,i)%UfeS = 0.0_dp
            evops(i,i)%Ugg = 1.0_dp

            if (.not.resources_have_next_block()) exit

            i = i + 1

            call resources_next_block()
            call resources_set_all_pointers(RESOURCES_EXCITON_REP)


        end do



    end subroutine allocate_evolution_operators



    !
    ! Calculates evolution operators
    !
    subroutine prepare_evolution_operators

        complex(dpc), dimension(:,:), pointer :: lUeg, fUeg
        integer :: i,k,l,t, kk
        real(dp) :: tt
        real(dp), dimension(:), pointer :: dr
        real(dp), dimension(:,:), pointer :: rr
        real(dp), dimension(:,:), allocatable :: KKr,Kexp,SS,S1
        integer :: N2
		character(len=64) :: cbuff

        complex(dpc), dimension(:,:,:), pointer :: gg_2
        complex(dpc), dimension(:,:,:), pointer :: gg_21


        ! ***********************************************************

		ggC = 1.0d0/100.0d0


        call resources_rewind_blocks()
        call resources_set_all_pointers(RESOURCES_EXCITON_REP)
        call rewind_evolution_operators()

        !call new_evolution_operators()

        write(cbuff,'(a,i3)') "Preparing evolution operators"
        call print_log_message(trim(cbuff),6)


        do

            ALLOCATE(lUeg,(N1,Nt(1)))
            ALLOCATE(fUeg,(N1,Nt(1)))

            dr => current_e_block%dr0

            ! loop over time
            do t = 1, Nt(1)

                tt = (t-1)*gt(1)*dt

                do k = 1, N1

                    kk = k

                    lUeg(k,t) = exp(-gg(kk,t)-(0.0_dp,1.0_dp)*(en(k)-rwa)*tt - dr(k)*(tt + (exp(-ggC*tt)-1.0d0)/ggC))
                    fUeg(k,t) = exp(-conjg(gg(kk,t))-(0.0_dp,1.0_dp)*(en(k)-rwa-2.0_dp*ll(k))*tt - dr(k)*tt)

                end do


            end do

            current_Ueg%Ueg = lUeg
            current_Ueg%feg = fUeg

			DEALLOCATE(lUeg)
			DEALLOCATE(fUeg)

			! relaxation constants
			rr => iblocks(1,1)%eblock%rr
			ALLOCATE(KKr,(N1,N1))
			ALLOCATE(Kexp,(N1,N1))
			ALLOCATE(SS,(N1,N1))
			ALLOCATE(S1,(N1,N1))

			call spec(rr,SS,KKr)

			call inv(SS,S1)

			! exciton block

			! evops(1,1)%Uee(k,l,k,l,t)
			! loop over time
            do t = 1, Nt(2)

                tt = (t-1)*gt(2)*dt

				Kexp = 0.0_dp
				if (ggC > 0.0d0) then
					do i = 1,N1
						Kexp(i,i) = exp(-KKr(i,i)*(tt + (exp(-ggC*tt)-1.0d0)/ggC))
					end do
				else
					do i = 1,N1
						Kexp(i,i) = exp(-KKr(i,i)*tt)
					end do
				end if

				Kexp = matmul(SS,matmul(Kexp,S1))

                do k = 1, N1
    			do l = 1, N1

					if (k /= l) then

						! dephasing
						evops(1,1)%Uee(k,l,k,l,t) = exp(-(0.0d0,1.0d0)*(en(k)-en(l))*tt - (KKr(l,l)+KKr(k,k))*tt )

					end if

					! relaxation
					evops(1,1)%Uee(k,k,l,l,t) = Kexp(k,l)

				   ! ! depopulation rates
				   ! evops(1,1)%Uee(k,k,k,k,t) = Kexp(k,k)


				end do
				end do


			end do



			DEALLOCATE(KKr)
			DEALLOCATE(SS)
			DEALLOCATE(S1)


			! one- to two-exciton coherences

            if (use_twoexcitons) then

				gg_2 => iblocks(1,1)%eblock%gg_2
				gg_21 => iblocks(1,1)%eblock%gg_21

				N2 = N1*(N1-1)/2
				! loop over time
            	do t = 1, Nt(1)

                	tt = (t-1)*gt(1)*dt

                	do k = 1, N2
    				do l = 1, N1

						! this is in fact Uef !!!!!
    		!*			evops(1,1)%Ufe(k,l,k,l,t) = exp(-conjg(gg_2(k,k,t)) - (gg(l,t))  + 2*real(gg_21(k,l,t)) &
    		!*			+(0.0d0,1.0d0)*( iblocks(1,1)%eblock%en_2(k) - iblocks(1,1)%eblock%en(l) - rwa)*tt - 2*dr(l)*tt)

    					evops(1,1)%UfeS(k,l,t) = exp(-conjg(gg_2(k,k,t)) - (gg(l,t))  + 2*real(gg_21(k,l,t)) &
    					+(0.0d0,1.0d0)*( iblocks(1,1)%eblock%en_2(k) - iblocks(1,1)%eblock%en(l) - rwa)*tt - 2*dr(l)*(tt + (exp(-ggC*tt)-1.0d0)/ggC))

            	    end do
	            	end do


                end do

            end if


            if (.not.resources_have_next_block()) exit

            call resources_next_block()
            call resources_set_all_pointers(RESOURCES_EXCITON_REP)
            call new_evolution_operators()


        end do

        ! ***************************************************************************

if (.false.) then

        ! new scheme of storing evolution operator
        call resources_rewind_blocks()
        call resources_set_all_pointers(RESOURCES_EXCITON_REP)
        !allocate(dr(N1))

        i = 1

        do

            N2 = N1*(N1-1)/2
            evops(i,i)%Ueg = 0.0_dp
            evops(i,i)%Ugg = 1.0_dp

            if (.false.) then

            dr => current_e_block%dr0

            ! loop over time
            do t = 1, Nt(1)

                tt = (t-1)*gt(1)*dt

                do k = 1, N1

                    !kk = k + (k-1)*N1
                    kk = k

                    evops(i,i)%Ueg(k,1,k,1,t) = exp(-gg(kk,t)-(0.0_dp,1.0_dp)*(en(k)-rwa)*tt - dr0(k)*(tt + (exp(-ggC*tt)-1.0d0)/ggC))
                    !fUeg(k,t) = exp(-conjg(gg(kk,t))-(0.0_dp,1.0_dp)*(en(k)-rwa-2.0_dp*ll(k))*tt - dr0(k)*tt)

                end do

            end do
            end if


            if (.not.resources_have_next_block()) exit
            i = i + 1

            call resources_next_block()
            call resources_set_all_pointers(RESOURCES_EXCITON_REP)


        end do

end if

    end subroutine prepare_evolution_operators

    !
    !
    !
    subroutine rewind_evolution_operators()
        current_Ueg => storage
    end subroutine rewind_evolution_operators

    !
    !
    !
    subroutine new_evolution_operators()
        if (current_Ueg%id == 0) then
            current_Ueg%id = 1
            current_Ueg%next => null()
        else
            allocate(current_Ueg%next)
            current_Ueg%next%id = current_Ueg%id + 1
            current_Ueg => current_Ueg%next
            current_Ueg%next=>null()
        end if
    end subroutine new_evolution_operators

    !
    !
    !
    subroutine next_evolution_operators()
        current_Ueg => current_Ueg%next
    end subroutine next_evolution_operators

    !
    ! Cleans all stuff
    !
    subroutine clean_resources_tdpt3()
        type(e_ops_storage), pointer              :: pp

        call rewind_evolution_operators()

        do


            if (associated(current_Ueg%Ueg)) then
                deallocate(current_Ueg%Ueg)
                !deallocate(current_Ueg%feg)
            end if

            pp => current_Ueg

            if (.not.associated(current_Ueg%next)) then
                exit
            end if

            call next_evolution_operators()

            deallocate(pp)


        end do

        deallocate(pp)

        !if (associated(storage)) deallocate(storage)
        !storage => null()

        ! initialize excion blocks
        !allocate(storage)
        !storage%id = 0
        !storage%next => null()
        !storage%Ueg  => null()

        !if (associated(Ueg)) then
        !    deallocate(Ueg)
        !end if
    end subroutine clean_resources_tdpt3


    subroutine super_unity(N,Nt,sI)
        integer :: N, Nt
        complex(dpc), dimension(:,:,:,:,:) :: sI

        integer :: i,j

        sI = 0.0_dp

        do i = 1, N
            do j = 1, N
                sI(i,j,i,j,1:Nt) = 1.0_dp
            end do
        end do


    end subroutine super_unity

end module resources_tdpt3
