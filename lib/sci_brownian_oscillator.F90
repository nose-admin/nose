!
! Stochastic Brownian oscilator
!
! Correlation functions and their up to second integral (g of t) 
!
! Author: Tomas Mancal, tomas.mancal@mff.cuni.cz
!
!
module sci_brownian_oscillator

	use std_types
	use numer_interp     ! interpolations, integrations

	implicit none
	
	real(dp), private     :: T
	real(dp), private     :: LL
	real(dp), private     :: l_r
	real(dp), private     :: hbeta
    real(dp), private     :: om, gamma
    complex(dpc), private :: phi
    real(dp), private     :: dzeta
	real(dp), private, dimension(:), allocatable :: nu
	! energy conversion
	real(dp), private     :: ConvE
	integer(i4b), private :: N_mats 
	logical, private :: initiated
    
    integer, private :: bo_tp
	
	private :: init_bo_N, init_bo_noarg
	
	public  :: init_brownian_oscillator 
	public  :: clean_brownian_oscillator
    public :: brownian_oscillator_set_type
    public :: brownian_oscillator_set_temperature
    public :: brownian_oscillator_set_frequency
    public :: brownian_oscillator_set_gamma
	
	interface init_brownian_oscillator
		module procedure init_bo_N, init_bo_noarg
	end interface 
    
    interface brownian_oscillator_set_type
        module procedure set_type
    end interface

    interface brownian_oscillator_set_temperature
        module procedure set_temperature
    end interface

    interface brownian_oscillator_set_frequency
        module procedure set_frequency
    end interface

    interface brownian_oscillator_set_correlation_time
        module procedure set_correlation_time
    end interface
    
    interface brownian_oscillator_set_gamma
        module procedure set_gamma
    end interface

    
    private :: set_type
    private :: set_temperature
    private :: set_frequency
    private :: set_correlation_time
    private :: set_gamma
	
contains

	!
	!
	!
	subroutine init_bo_N(Nin)
		integer(i4b), intent(in) :: Nin
    	if (initiated) stop
		ConvE  = 2.0_dp*acos(-1.0_dp)*2.99792458d-5
		N_mats = Nin
		allocate(nu(N_mats))
		T      = 100.0_dp
		LL     = 1.0_dp/1000.0_dp
		hbeta  = (0.6582120_dp/8.617385d-5)/T  ! hbar/kT
		call matsubara()
		initiated = .true.
        bo_tp = 0
	end subroutine init_bo_N

	!
	!
	!
	subroutine init_bo_noarg()
          integer(i4b) :: N
          N = 100
          call init_bo_N(N)	
	end subroutine init_bo_noarg

	!
	!
	!
	subroutine clean_brownian_oscillator()
		deallocate(nu)
	end subroutine clean_brownian_oscillator


    subroutine brownian_oscillator_update()
        select case (bo_tp)
        case (0)    
        case (1)
          !gamma = (om**2)/LL
          dzeta = sqrt(om**2-(gamma/2.0_dp)**2)
          phi = gamma/2.0_dp + (0.0_dp,1.0_dp)*dzeta
        case default
            print *, "Unknown type"
        end select        
    end subroutine brownian_oscillator_update

    subroutine set_type(type)
        character(len=*), intent(in) :: type
        if (.not.initiated) stop
        
        if (type == "GENERAL") then
            bo_tp = 1
        else if (type == "OVERDAMPED") then
            bo_tp = 0
        else 
            print *, "Unknown BO type"
            stop    
        end if        
        
    end subroutine set_type

	!
	! sets the temperature in K
	!
	subroutine set_temperature(TT) 
		real(dp), intent(in) :: TT
		if (.not.initiated) stop
		T = TT
		if (T /= 0.0_dp) then
			hbeta  = (0.6582120_dp/8.617385d-5)/T  ! hbar/kT
		else 
			hbeta = 0.0_dp
		end if
		call matsubara()
	end subroutine set_temperature

    subroutine set_frequency(omi)
        real(dp), intent(in) :: omi
        if (.not. initiated) stop
        om = omi
    end subroutine set_frequency

	!
	! sets correlation time in fs
	!
	subroutine set_correlation_time(ttauc)
		real(dp), intent(in) :: ttauc
		if (.not.initiated) stop
		LL =  1.0_dp/ttauc
	end subroutine set_correlation_time
	
    subroutine set_gamma(gg)
        real(dp), intent(in) :: gg
        gamma = gg        
    end subroutine set_gamma
    
	!
	! calculates Matsubara frequencies
	!
	subroutine matsubara()
		integer(i4b) :: i
		do i = 1, N_mats
			nu(i) = 2.0_dp*PI_D*real(i,dp)/hbeta
		end do
	end subroutine matsubara
	

!
! BO CORRELATION FUNCTIONS IN TIME DOMAIN
!	
	
	!
	! C''(t) for reorganization energy equal 1 
	!
	function Cdd1_t(t) result (ret)
		real(dp), intent(in) :: t
		real(dp)             :: ret
	
        select case (bo_tp)
        case (0)    
          ret = -LL*exp(-LL*t)
        case (1)
          ret = -exp(-gamma*t/2.0_dp)*sin(dzeta*t)/(2.0_dp*dzeta)
        case default
            print *, "Unknown type"
        end select
        
            
	end function Cdd1_t
	
	!
	! C'(omega) for reorganization energy equal 1
	!
	function Cd1_t(t) result (ret)
		real(dp), intent(in) :: t
		real(dp)             :: ret
        
        select case (bo_tp)
        case (0)    
		  ret = (LL/tanh(hbeta*LL/2.0_dp))*exp(-LL*t) + &
			4.0_dp*LL*msum_t(t)/hbeta
        case (1)
          ret = ((1.0_dp/ctanh((0.0_dp,1.0_dp)*conjg(phi)*hbeta/2.0_dp))*exp(-conjg(phi)*t) - &
                 (1.0_dp/ctanh((0.0_dp,1.0_dp)*      phi* hbeta/2.0_dp))*exp(-       phi*t))/(4.0_dp*dzeta) - &
                 (2.0_dp*gamma/hbeta)*msum_gen_t(t)
        case default
            print *, "Unknown type"
        end select

!        print *, t, (0.0_dp,1.0_dp)*conjg(phi)*hbeta/2.0_dp, ctanh((0.0_dp,1.0_dp)*conjg(phi)*hbeta/2.0_dp) , &
!                                                             ctanh((0.0_dp,1.0_dp)*phi*hbeta/2.0_dp)
!                (1.0_dp/ctanh((0.0_dp,1.0_dp)*conjg(phi)*hbeta/2.0_dp)*exp(-conjg(phi)*t) - &
!                 1.0_dp/ctanh((0.0_dp,1.0_dp)*phi*hbeta/2.0_dp)*exp(-phi*t))/(4.0_dp*dzeta) - &
!                 (2.0_dp*gamma/hbeta)*msum_gen_t(t)


	end function Cd1_t

	!
	! Matsubara term in time domain for 
	! reorganization energy equal 1 
	!
	function msum_t(tt) result (ret)
		real(dp) :: tt
		real(dp) :: ret 
		integer(i4b) :: i
		ret = 0.0_dp
		if (T == 0.0_dp) return
		do i = 1, N_mats
			ret = ret + (nu(i)*exp(-nu(i)*tt))/  &
				(nu(i)**2 - LL**2) 
		end do
	end function msum_t
    
    function msum_gen_t(tt) result (ret)
        real(dp) :: tt
        real(dp) :: ret 
        integer(i4b) :: i
        ret = 0.0_dp
        if (T == 0.0_dp) return
        do i = 1, N_mats
            ret = ret + (nu(i)*exp(-nu(i)*tt))/  &
                ((om**2+nu(i)**2) - (gamma*nu(i))**2) 
        end do
    end function msum_gen_t
	

!
! NUMERICAL INTEGRATION OF BO CORRELATION FUNCTIONS
!

	subroutine Hdd1_t(hd,tt,omega) 
		complex(dpc), dimension(:), intent(out) :: hd
		real(dp), dimension(:), intent(in)      :: tt
		real(dp), intent(in)                    :: omega
	
		! local
		integer(i4b) :: sh,st, Nt, i
		complex(dpc), dimension(size(tt),1) :: pyr,pyi
		real(dp), dimension(size(tt)) :: yy 
		real(dp), dimension(1)              :: w
	
		sh = size(hd)
		st = size(tt)
		Nt = st
		
!		print *, st, sh
		
		if (sh /= st) stop
		
		w(1) = omega
		do i = 1, Nt
			yy(i) = Cdd1_t(tt(i))
		end do

        call primitive(tt,yy,w,pyr)

		hd = pyr(:,1)
	
	end subroutine Hdd1_t

    subroutine Gdd1_t(hd,tt,omega) 
        complex(dpc), dimension(:), intent(out) :: hd
        real(dp), dimension(:), intent(in)      :: tt
        real(dp), intent(in)                    :: omega
    
        ! local
        integer(i4b) :: sh,st, Nt, i
        complex(dpc), dimension(size(tt),1) :: pyr,pyi
        real(dp), dimension(size(tt)) :: yy 
        real(dp), dimension(1)              :: w
    
        sh = size(hd)
        st = size(tt)
        Nt = st
        
!       print *, st, sh
        
        if (sh /= st) stop
        
        w(1) = omega
        do i = 1, Nt
            yy(i) = Cdd1_t(tt(i))
        end do

        ! first integration
        call primitive(tt,yy,w,pyi)

        yy = real(pyi(:,1))
        
        ! second integration - real part
        call primitive(tt,yy,w,pyr)
                
        yy = aimag(pyi(:,1))

        ! second integration - imaginary part
        call primitive(tt,yy,w,pyi)

        hd = pyr(:,1) + (0.0_dp,1.0_dp)*pyi(:,1)
    
    end subroutine Gdd1_t
	

	subroutine Hd1_t(hd,tt,omega) 
		complex(dpc), dimension(:), intent(out) :: hd
		real(dp), dimension(:), intent(in)      :: tt
		real(dp), intent(in)                    :: omega
	
		! local
		integer(i4b) :: sh,st, Nt, i
		complex(dpc), dimension(size(tt),1) :: pyr,pyi
		real(dp), dimension(size(tt)) :: yy 
		real(dp), dimension(1)              :: w
	
		sh = size(hd)
		st = size(tt)
		Nt = st
		
!		print *, st, sh
		
		if (sh /= st) stop
		
		w(1) = omega
		do i = 1, Nt
			yy(i) = Cd1_t(tt(i))
		end do

        call primitive(tt,yy,w,pyr)

        hd = pyr(:,1)
	
	end subroutine Hd1_t


    subroutine Gd1_t(hd,tt,omega) 
        complex(dpc), dimension(:), intent(out) :: hd
        real(dp), dimension(:), intent(in)      :: tt
        real(dp), intent(in)                    :: omega
    
        ! local
        integer(i4b) :: sh,st, Nt, i
        complex(dpc), dimension(size(tt),1) :: pyr,pyi
        real(dp), dimension(size(tt)) :: yy 
        real(dp), dimension(1)              :: w
    
        sh = size(hd)
        st = size(tt)
        Nt = st
        
!        print *, st, sh
        
        if (sh /= st) stop
        
        w(1) = omega
        do i = 1, Nt
            yy(i) = Cd1_t(tt(i))
        end do

        ! first integration
        call primitive(tt,yy,w,pyi)

        yy = real(pyi(:,1))
        
        ! second integration - real part
        call primitive(tt,yy,w,pyr)
                
        yy = aimag(pyi(:,1))

        ! second integration - imaginary part
        call primitive(tt,yy,w,pyi)

        hd = pyr(:,1) + (0.0_dp,1.0_dp)*pyi(:,1)
    
    end subroutine Gd1_t

	!
	! Numerical Half sided FT of C''(t)
	!
	function Hdd1_inf(omega) result (ret)
		real(dp), intent(in)                    :: omega
		complex(dpc) :: ret 
		
		! local
		integer(i4b), parameter :: Nt = 2000
		complex(dpc), dimension(Nt) :: hd
		real(dp), dimension(Nt)     :: tt
		integer(i4b) :: i
		complex(dpc), dimension(Nt,1) :: pyr,pyi
		real(dp), dimension(Nt) :: yy 
		real(dp), dimension(1)              :: w
			
		w(1) = omega
		do i = 1, Nt
			tt(i) = real(i-1,dp)
			yy(i) = Cdd1_t(tt(i))
		end do

        call primitive(tt,yy,w,pyr)

		ret = pyr(Nt,1)
	
	end function Hdd1_inf
	
	function Hd1_inf(omega) result (ret)
		real(dp), intent(in)                    :: omega
		complex(dpc) :: ret
	 
	 	integer(i4b), parameter :: Nt = 2000
		complex(dpc), dimension(Nt) :: hd
		real(dp), dimension(Nt)     :: tt
		integer(i4b) :: i
		complex(dpc), dimension(size(tt),1) :: pyr,pyi
		real(dp), dimension(size(tt)) :: yy 
		real(dp), dimension(1)              :: w
	
		w(1) = omega
		do i = 1, Nt
			tt(i) = real(i-1,dp)
			yy(i) = Cd1_t(tt(i))
		end do

        call primitive(tt,yy,w,pyr)

		ret = pyr(Nt,1)
	
	end function Hd1_inf
		
!
! HALF SIDE FOURIER TRANSFORMED BO CORRELATION FUNCTIONS	
!	
	
	!
	! C''(omega) for reorganization energy equal 1 
	!
	function Cdd1(omega) result (ret)
		real(dp), intent(in) :: omega
		complex(dpc) :: ret
		ret = -(omega*LL + (0.0_dp,1.0_dp)*LL**2)/(omega**2 + LL**2)
	end function Cdd1
	
	!
	! C'(omega) for reorganization energy equal 1
	!
	function Cd1(omega) result (ret)
		real(dp) :: omega
		complex(dp) :: ret
		if (T == 0.0_dp) then
			ret = 0.0_dp
			return
		end if 
		ret = (LL**2 - (0.0_dp,1.0_dp)*omega*LL)/ &
			((omega**2 + LL**2)*tanh(hbeta*LL/2.0_dp)) &
			+ 4.0_dp*LL*msum(omega)/hbeta
	end function Cd1

	!
	! Half side Fourier transformed Matsubara term for 
	! reorganization energy equal 1 
	!
	function msum(omega) result (ret)
		real(dp) :: omega
		complex(dpc) :: ret 
		integer(i4b) :: i
		ret = 0.0_dp
		if (T == 0.0_dp) return
		do i = 1, N_mats
			ret = ret + (nu(i)**2-(0.0_dp,1.0_dp)*omega*nu(i))/  &
				((nu(i)**2 - LL**2)*(nu(i)**2 + omega**2)) 
		end do
	end function msum
    
    function ctanh(x) result (z)
        complex(dpc), intent(in) :: x
        complex(dpc)             :: z
        
        z = (exp(x)-exp(-x))/(exp(x)+exp(-x))     
        
    end function ctanh 
	
	
end module sci_brownian_oscillator
