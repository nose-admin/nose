!
!  Module handling calculation of the 3rd order signal using response functions
!
!
module signal3

    use std_types
#include <util_allocation.h>
    use util_allocation
    use response3

    complex(dpc), dimension(:,:), allocatable, public :: SR0, SNR0
    complex(dpc), dimension(:,:), allocatable, public :: ESAN,ESAR


contains


    subroutine init_signal3()

        integer :: Nt1, Nt2

		nonsecular = .false.

        call init_response3()
        NFFT_basic = 1024

        padfac = 0
        NFFT = (2**padfac)*NFFT_basic
        Nt1 = size(R2g0,1)
        Nt2 = size(R2g0,2)


        ALLOCATE(SR0,(NFFT,NFFT))
        ALLOCATE(SNR0,(NFFT,NFFT))
        ALLOCATE(ESAR,(NFFT,NFFT))
        ALLOCATE(ESAN,(NFFT,NFFT))

        SR0 = 0.0_dp
        SR0(1:Nt1,1:Nt2) = SR0(1:Nt1,1:Nt2) + R2g0
        SR0(1:Nt1,1:Nt2) = SR0(1:Nt1,1:Nt2) + R3g0
        !SR0(1:Nt1,1:Nt2) = SR0(1:Nt1,1:Nt2) - conjg(R1f0)

        SNR0 = 0.0_dp
        SNR0(1:Nt1,1:Nt2) = SNR0(1:Nt1,1:Nt2) + R1g0
        SNR0(1:Nt1,1:Nt2) = SNR0(1:Nt1,1:Nt2) + R4g0
        !SNR0(1:Nt1,1:Nt2) = SNR0(1:Nt1,1:Nt2) - conjg(R2f0)

		ESAR = 0.0_dp
		ESAR(1:Nt1,1:Nt2) = ESAR(1:Nt1,1:Nt2) - conjg(R1f0)
		ESAN = 0.0_dp
		ESAN(1:Nt1,1:Nt2) = ESAN(1:Nt1,1:Nt2) - conjg(R2f0)

    end subroutine init_signal3


    subroutine clean_signal3()
        DEALLOCATE(SR0)
        DEALLOCATE(SNR0)
        DEALLOCATE(ESAN)
        DEALLOCATE(ESAR)
        call clean_response3()
    end subroutine clean_signal3


end module signal3
