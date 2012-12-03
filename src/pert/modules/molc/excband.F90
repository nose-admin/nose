module excband
    
    use std_types
    
    
    
contains


    subroutine excband_structure()
        integer :: Nx
        integer :: i

        real(dp), dimension(:), allocatable :: ee
        complex(dp), dimension(:), allocatable :: dd

        write(*,*) 'Linear chain calculation'
        
        Nx = 100
        allocate(ee(Nx),dd(Nx))
        
        ! energy
        do i = 1, Nx
            
            ee(i) = 0.0_dp + 2*J*cos((TwoPi_D/Nx)*i)
            
        end do
        
        open(unit=11,file="energy.dat")
        do i = 1, Nx
            write(11,*) i, ee(i)    
        end do
        close(11)
        
        ! dipole moment
        do i = 1, Nx
            dd(i) = 0.0_dp
            do j = 1, Nx
                
                dd(i) = dd(i) + cos(j*TwoPi_D/(2*Nx-3))*exp(-(0.0_dp,1.0_dp)*(TwoPi_D/Nx)*i*j)
                
            end do 
            
        end do
        
        open(unit=11,file="dsquare.dat")
        do i = 1, Nx
            write(11,*) i, abs(dd(i))**2    
        end do
        close(11)
        
        ! absorption spectrum
        
        
    end subroutine excband_structure
    
    
    
    
    
end module excband
