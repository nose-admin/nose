module sci_liouville
   
    use std_types
    use numer_random
    
    implicit none
   
    real(dp), private, dimension(:), allocatable :: gamma
    real(dp), private, dimension(:), allocatable :: omega
    real(dp), private                            :: dt
    real(dp), private, dimension(:), allocatable :: om0
        
    
    complex(dpc), private, dimension(:), allocatable :: KK
    real(dp), private, dimension(:), allocatable :: a2,b2,ab,cc  ! moments of the Gaussian distribution
                                                                 ! and correlation
    integer  :: N_mods
   
    real(dp), private :: Temp ! temperature
    
    integer, private  :: par_count
    
    ! energy conversion factor cm-1 -> internal 
    !real(dp), parameter :: Energy_cm_to_internal = 2.0_dp*PI_D*2.99792458d-5
    !real(dp), parameter :: Energy_internal_to_cm = 5308.8415534785_dp
   
   
contains
    !
    !
    !
    subroutine init_liouville(N, T)
        integer, intent(in)  :: N  ! number of modes
        real(dp), intent(in) :: T
    
        allocate(gamma(N),omega(N),om0(N),KK(N),a2(N),b2(N),ab(N),cc(N))
        gamma = 0.0_dp
        omega = 0.0_dp
        om0   = 0.0_dp
        dt    = 1.0_dp
        Temp  = T
        par_count = 0
        N_mods = N
        
    end subroutine init_liouville  
    
    !
    ! Set all gammas
    !
    subroutine liouville_set_gamma(g)
        real(dp), dimension(:), intent(in) :: g
        
        if (size(g)==N_mods) then
            gamma = g
            par_count = par_count + 1
        else
            stop
        end if
        
        if (par_count >= 4) then
            call reset_moments()
        end if
        
    end subroutine liouville_set_gamma
    
    subroutine liouville_get_gamma(g)
        real(dp), dimension(:), intent(out) :: g
        g = gamma
    end subroutine liouville_get_gamma


    !
    ! Set all omegas
    !
    subroutine liouville_set_omega(g)
        real(dp), dimension(:), intent(in) :: g
        
        if (size(g)==N_mods) then
            omega = g
            par_count = par_count + 1
        else
            stop
        end if
        
        if (par_count >= 4) then
            call reset_moments()
        end if
        
    end subroutine liouville_set_omega
    
    subroutine liouville_get_omega(g)
        real(dp), dimension(:), intent(out) :: g
        g = omega
    end subroutine liouville_get_omega
    
    !
    ! Set all om0s
    !
    subroutine liouville_set_om0(g)
        real(dp), dimension(:), intent(in) :: g
        
        if (size(g)==N_mods) then
            om0 = g
            par_count = par_count + 1
        else
            stop
        end if
               
        if (par_count >= 4) then
            call reset_moments()
        end if
        
    end subroutine liouville_set_om0
    
    subroutine liouville_get_om0(g)
        real(dp), dimension(:), intent(out) :: g
        g = om0
    end subroutine liouville_get_om0
    
    !
    ! Set dt
    !
    subroutine liouville_set_dt(ddt)
        real(dp),intent(in) :: ddt
        
        dt = ddt
        par_count = par_count + 1
        
        if (par_count >= 4) then
            call reset_moments()
        end if
        
    end subroutine liouville_set_dt
    
    subroutine liouville_get_dt(g)
        real(dp), intent(out) :: g
        g = dt
    end subroutine liouville_get_dt

    
    
    
    !
    ! Rest the moments of the distribution
    !
    subroutine reset_moments()
        integer  :: i
        complex(dpc) :: pom
        real(dp) :: rpom
        real(dp) :: X, Y, Z
        real(dp) :: kb
        kb = Energy_cm_to_internal*kB_cmK
        do i = 1, N_mods
            rpom = omega(i)**2 - (gamma(i)/2)**2 
            !print *, rpom
            !if (rpom >= 0.0_dp) then
                KK(i) = gamma(i)/2.0_dp - (0.0_dp,1.0_dp)*sqrt(rpom)/2.0_dp
            !else
            !    KK(i) = gamma(i)/2.0_dp + sqrt(-rpom)/2.0_dp
            !end if
            !print *, real(KK)
            pom  = 2.0_dp*gamma(i)*kb*Temp*(1.0_dp-exp(-2.0_dp*KK(i)*dt))/(2.0_dp*KK(i)) 
            X = real(pom)
            Y = aimag(pom)
            Z = 2.0_dp*gamma(i)*kb*Temp*(1.0_dp-exp(-2.0_dp*real(KK(i))*dt))/(2.0_dp*real(KK(i)))
            a2(i) = X + (X+Z)/2.0_dp
            b2(i) = (X+Z)/2.0_dp
            ab(i) = Y/2.0_dp
            cc(i) = ab(i)/sqrt(a2(i)*b2(i))
            !print *, "Correlation: ", cc(i) 
        end do
    end subroutine reset_moments
    
    subroutine clean_liouville()
        deallocate(gamma,omega,om0,KK,a2,b2,ab)        
    end subroutine clean_liouville
    
    
    subroutine liouville_make_step(Q)
        complex(dpc), dimension(:), intent(inout) :: Q
        integer :: i
        if (size(Q)==N_mods) then
            do i = 1, N_mods
                Q(i) = exp(-KK(i)*dt)*Q(i) + FF(i)  + (1.0_dp-exp(-KK(i)*dt))*om0(i)/KK(i) 
            end do
        else
            stop
        end if
    end subroutine liouville_make_step
    
    subroutine liouville_get_coordinates(Q,qq)
        complex(dpc), dimension(:), intent(in) :: Q
        real(dp), dimension(:), intent(out)    :: qq
        integer :: i
        
        if (size(Q)==N_mods) then
            do i = 1, N_mods
                qq(i) = aimag(Q(i))/aimag(gamma(i)-KK(i))
            end do
        else
            stop
        end if       
    end subroutine liouville_get_coordinates 

    subroutine liouville_get_impulses(Q,qq)
        complex(dpc), dimension(:), intent(in) :: Q
        real(dp), dimension(:), intent(out)    :: qq
        complex(dpc) :: xi
        integer :: i
        
        if (size(Q)==N_mods) then
            do i = 1, N_mods
                xi = gamma(i)-KK(i)
                qq(i) = (xi*conjg(Q(i))-conjg(xi)*Q(i))/((0.0_dp,2.0_dp)*aimag(xi))
            end do
        else
            stop
        end if       
    end subroutine liouville_get_impulses


    
    subroutine liouville_get_xi(xi)
        complex(dpc), dimension(:), intent(out) :: xi
        xi(:) = gamma(:) - KK(:)    
    end subroutine liouville_get_xi
    
    !
    ! Random force
    !
    function FF(i) result (f)
        integer, intent(in) :: i
        complex(dpc) :: f
        real(dp) :: a,b
        call biNormal(a,b,i)
        f =  a + (0.0_dp,1.0_dp)*b
        !f = 0.0_dp
    end function FF
    
    
    !
    !
    !
    subroutine biNormal(a,b,i) 
        real(dp), intent(out)  :: a,b
        integer,  intent(in)   :: i
        real(dp), dimension(2) :: xx
        
        call random_Normal(xx)
        
        a = sqrt(a2(i))*xx(1)
        b = cc(i)*sqrt(b2(i))*xx(1) + sqrt(1-cc(i)**2)*sqrt(b2(i))*xx(2)
        
    end subroutine biNormal
    
end module sci_liouville
