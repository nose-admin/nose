!
!  Module handling calculation of the 3rd order signal using response functions
!
!
module signal3

    use std_types
#include <util_allocation.h>
    use util_allocation
    use response3

    complex(dpc), dimension(:,:), allocatable, public :: SR0, SNR0, SSR1, SSR2
    complex(dpc), dimension(:,:), allocatable, public :: ESAN,ESAR


contains


    subroutine init_signal3()

        integer :: Nt1, Nt2

		write(*,*) modname
		if(modname == "QME") then
			if(sum(current_s_block%QHO_lvls) == 0) then
				nonsecular = .true.
				!use_Uee_index = INT(gt(2)/gt(1)) + 1
			else
				!use_Uee_index = 2
				nonsecular = .false.
				call init_orfact()
			end if
		elseif(modname == 'TDPT-3') then
			nonsecular = .false.
		else
			nonsecular = .true.
			call print_warning_message("It was not set whether to use secular EVOPS for this module, using non-secular.",-1)
		end if


        call init_response3()
        NFFT_basic = 1024

        padfac = 0
        NFFT = (2**padfac)*NFFT_basic
        Nt1 = size(R2g0,1)
        Nt2 = size(R2g0,2)
        write(*,*) 'padfac=',padfac

        ALLOCATE(SR0,(Nt1,Nt2))
        ALLOCATE(SNR0,(Nt1,Nt2))
        ALLOCATE(SSR1,(Nt1,Nt2))
        ALLOCATE(SSR2,(Nt1,Nt2))
        ALLOCATE(ESAR,(Nt1,Nt2))
        ALLOCATE(ESAN,(Nt1,Nt2))

        SR0 = 0.0_dp
        SR0(1:Nt1,1:Nt2) = SR0(1:Nt1,1:Nt2) + R2g0
        SR0(1:Nt1,1:Nt2) = SR0(1:Nt1,1:Nt2) + R3g0

        SSR1 = 0.0_dp
        SSR1(1:Nt1,1:Nt2) = SSR1(1:Nt1,1:Nt2) + R1g0

        SSR2 = 0.0_dp
        SSR2(1:Nt1,1:Nt2) = SSR2(1:Nt1,1:Nt2) + R2g0

        SNR0 = 0.0_dp
        SNR0(1:Nt1,1:Nt2) = SNR0(1:Nt1,1:Nt2) + R1g0
        SNR0(1:Nt1,1:Nt2) = SNR0(1:Nt1,1:Nt2) + R4g0

		ESAR = 0.0_dp
		ESAR(1:Nt1,1:Nt2) = ESAR(1:Nt1,1:Nt2) - conjg(R1f0)
		ESAN = 0.0_dp
		ESAN(1:Nt1,1:Nt2) = ESAN(1:Nt1,1:Nt2) - conjg(R2f0)

    end subroutine init_signal3


    subroutine clean_signal3()
        DEALLOCATE(SR0)
        DEALLOCATE(SNR0)
        DEALLOCATE(SSR1)
        DEALLOCATE(SSR2)
        DEALLOCATE(ESAN)
        DEALLOCATE(ESAR)
        call clean_response3()
    end subroutine clean_signal3


end module signal3
