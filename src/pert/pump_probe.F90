!
! Module containing subroutines for pump-probe spectroscopy
!
module pump_probe

    use std_types
    use signal3
    use numer_fft

contains


    subroutine init_pump_probe()

        integer :: i,j

!        Nt1 = Nt(1)
!
!
        call init_signal3()
        
        ! here SR0 and SNR0 contain the polarization as a function of t1 and t3.
        ! Pump-probe spectrum is calculated from the polarization at P(t3,t1=0)
        
        
        

        print *, 'PUMP-PROBE with SIGNAL'

        !stop

    end subroutine init_pump_probe

    
end module pump_probe