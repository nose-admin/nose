#include "util_allocation.h"
module resources_qme

    use resources
    use sci_redfield

    implicit none

    type populations
        integer                                 :: N
        real(dp), dimension(:,:), pointer       :: P
    end type

    type grstate_coherences
        integer                                 :: N
        complex(dpc), dimension(:,:), pointer   :: C
    end type

    type remaining_coherences
        integer                                   :: N
        complex(dpc), dimension(:,:,:), pointer   :: RC
    end type

    real(dp), dimension(:,:), allocatable :: K   ! matrix of pop. transfer rates
    complex(dpc), dimension(:,:,:), allocatable :: cpl ! matrix of couplings between sites in a block
    complex(dpc), dimension(:,:,:,:,:), allocatable :: RelMa ! the ralaxation matrix


    complex, dimension(:), allocatable :: qme_polar_1
    real, dimension(:), allocatable    :: qme_spect_abs

    complex, dimension(:,:), allocatable :: qme_2D_ftpe

    logical :: use_module_nakajima_zwanzig, tau_projector_normalization_for_others, &
    			qme_module_already_initiated = .false.

    type(populations)        :: pops
    type(grstate_coherences) :: gcohs
    type(remaining_coherences) :: rcohs


	integer, dimension(10) :: external_data
    character(len=64) :: method
    character          :: submethod1
    character          :: submethod2

!PGF95    private :: populations, grstate_coherences, remaining_coherences

contains

    !
    ! Initialization of the QME specific variables
    !
    subroutine init_resources_qme()

        integer :: N

 !       real(dp),dimension(:,:), allocatable :: HH

        N = N_ExOrder_Added(nr_blocks,1)

        pops%N  = N
        gcohs%N = N
        rcohs%N = N
        allocate(pops%P(Nt(1),N))
        allocate(gcohs%C(Nt(1),N))
        allocate(rcohs%RC(Nt(1),N,N))

!        allocate(HH(N,N))

!        HH = 0.0_dp

        call prepare_evolution_operators()

        external_data       = 0
        
        if (read_external_evops) then
            external_data(1:3) = 1
        end if

    end subroutine init_resources_qme


    !
    ! Cleaning the QME specific variables
    !
    subroutine clean_resources_qme()
		integer :: i

        deallocate(pops%P,gcohs%C,rcohs%RC)
        if (allocated(qme_polar_1)) then
            deallocate(qme_polar_1)
        end if
        if (allocated(qme_spect_abs)) then
            deallocate(qme_spect_abs)
        end if
        if (allocated(K)) then
            deallocate(K)
        end if
		if (allocated(cpl)) then
            deallocate(cpl)
        end if
        if (allocated(RelMa)) then
            deallocate(RelMa)
        end if
        ! undoing the prepare_evolution_operators()
                i = 1
        do
        	if (associated(evops(i,i)%Ueg)) then
       		DEALLOCATE(evops(i,i)%Ueg)
        	end if
       		if (associated(evops(i,i)%Uee)) then
       		DEALLOCATE(evops(i,i)%Uee)
        	end if
        	if (associated(evops(i,i)%Ufe)) then
       		DEALLOCATE(evops(i,i)%Ufe)
        	end if
        	if (associated(evops(i,i)%Ugg)) then
       		DEALLOCATE(evops(i,i)%Ugg)
        	end if

            if (.not.resources_have_next_block()) exit
           i = i + 1

        end do

    end subroutine clean_resources_qme


    subroutine prepare_evolution_operators()

        real(dp), dimension(:), pointer :: dr
        integer :: i,t,k,kk, N2 ! this one to be moved
        real(dp) :: tt

       ! new scheme of storing evolution operator
        call resources_rewind_blocks()
        call resources_set_all_pointers(RESOURCES_EXCITON_REP)
        !allocate(dr(N1))

        i = 1
        do

            N2 = N1*(N1-1)/2
            ALLOCATE(evops(i,i)%Ueg,(N1, 1,N1, 1,Nt(1)))
            ALLOCATE(evops(i,i)%Uee,(N1,N1,N1,N1,Nt(1)))
            ALLOCATE(evops(i,i)%Ufe,(N2,N1,N2,N1,Nt(1)))
            ALLOCATE(evops(i,i)%Ugg,(1,1,1,1,Nt(1)))
            evops(i,i)%Ueg = 0.0_dp
            evops(i,i)%Uee = 0.0_dp
            evops(i,i)%Ufe = 0.0_dp
            evops(i,i)%Ugg = 1.0_dp

            dr => current_e_block%dr0

			if (.false.) then
            ! loop over time
            do t = 1, Nt(1)

                tt = (t-1)*gt(1)*dt

                do k = 1, N1

                    !kk = k + (k-1)*N1
                    kk = k

                    evops(i,i)%Ueg(k,1,k,1,t) = exp(-gg(kk,t)-(0.0_dp,1.0_dp)*(en(k)-rwa)*tt - dr0(k)*tt)
                    !fUeg(k,t) = exp(-conjg(gg(kk,t))-(0.0_dp,1.0_dp)*(en(k)-rwa-2.0_dp*ll(k))*tt - dr0(k)*tt)

                end do

            end do
			end if


            if (.not.resources_have_next_block()) exit
            i = i + 1

            call resources_next_block()
            call resources_set_all_pointers(RESOURCES_EXCITON_REP)



        end do

    end subroutine prepare_evolution_operators


end module resources_qme
