!
!  Foerster energy transfer rates
!
!  Tomas Mancal, tomas.mancal@mff.cuni.cz
!
!
!  TODO: Pridat set a get methods a umoznit prepocitani vysledku se zmenami
!
!
module sci_foerster

	use numer_interp     ! interpolations, integrations
	               ! propagates: scilab, std_types
	
	implicit none
		
	public :: init_foerster, clean_foerster
	
	public :: foerster_rate, foerster_overlap
	public :: foerster_NA
	public :: foerster_ND

	! interfaces to private subroutines 
	interface init_foerster
		module procedure f_init_1
	end interface

	real(dp), dimension(:,:), allocatable :: foerster_rate, foerster_overlap
	integer(i4b)                          :: foerster_NA, foerster_ND

    !
    ! Private stuff
    !

	private :: f_init_1

contains

	!
	! Initialization for matrix
	!
	subroutine f_init_1(NA,ND,JJ,ggA,ggD,llA,llD,Nt,dt,Om,Temp)
		integer,                 intent(in)   :: NA
		integer,                 intent(in)   :: ND
		real(dp),     dimension(:,:), intent(in)   :: JJ
		complex(dp),  dimension(:,:), intent(in)   :: ggA
		complex(dp),  dimension(:,:), intent(in)   :: ggD
		real(dp),     dimension(:),   intent(in)   :: llA
		real(dp),     dimension(:),   intent(in)   :: llD
		integer,                        intent(in)   :: Nt
		real(dp),                     intent(in)   :: dt
		real(dp),     dimension(:,:), intent(in)   :: Om
		real(dp),                     intent(in)   :: Temp

        ! local
		real(dp)                                   :: integ
		real(dp),     dimension(Nt)                :: tt
		real(dp),     dimension(1)                 :: w
		complex(dpc), dimension(Nt)                :: yy,yy2
		complex(dpc), dimension(Nt,1)              :: pyr,pyi
		integer(i4b)                               :: Ntest
		integer(i4b)                               :: i,j
		integer(i4b)                               :: mu,nu
		real(dp)                                   :: wi, wc
		real(dp)                                   :: convET,convFT, bal
		real(dp)                                   :: OmIn
		logical                                    :: bal_test
		
		bal_test = .false.

   		ConvET = 2.0_dp*acos(-1.0_dp)*2.99792458d-5
    	ConvFT = 1.43877_dp/Temp
				
		! check if the rate matrix is allocated (deallocate if it is)
		if (allocated(foerster_rate)) then
			deallocate(foerster_rate)
			deallocate(foerster_overlap)
		end if
		
		! allocation
		allocate(foerster_rate(NA,ND),foerster_overlap(NA,ND))
		foerster_NA = NA
		foerster_ND = ND
		 
		! check the input
		if (size(JJ,1)  /= NA) stop "JJ 1"
		if (size(JJ,2)  /= ND) stop "JJ 2"
		if (size(Om,1)  /= NA) stop "Om 1"
		if (size(Om,2)  /= ND) stop "Om 2"
		if (size(ggA,2) /= NA) stop "ggA"
		if (size(ggD,2) /= ND) stop "ggB"
		
		! create time line
		do i = 1, Nt
			tt(i) = real(i-1,dp)*dt
		end do
		
		! clear rates
		foerster_rate=0.0_dp
		 
		! calculate rates
		do mu = 1, NA
			do nu = 1, ND

            	w(1) = -Om(mu,nu) !-Om(mu,nu)
            	
!            	print *, "ome DA: ",w(1)
            	
                do i = 1, Nt
                	yy(i) = exp(-ggA(i,mu)-ggD(i,nu) -  &
                              (0.0_dp,1.0_dp)*(llA(mu)+llD(nu))*tt(i))
                end do

                call primitive(tt,real(yy),w,pyr)
                call primitive(tt,aimag(yy),w,pyi)

                integ = 2.0_dp*real(pyr(Nt,1)+(0.0_dp,1.0_dp)*pyi(Nt,1))

                ! correction for numerical or other errors
                ! always check the foerster.log to see the significance of
                ! negative rates
                if (integ > 0.0_dp) then
                	foerster_rate(mu,nu) =  &
                		(abs(JJ(mu,nu))**2)*integ
                    foerster_overlap(mu,nu) = integ
                else
                	foerster_rate(mu,nu) = 0.0_dp
                    foerster_overlap(mu,nu) = 0.0_dp
                end if
								
			end do		
		end do
	
	
		!
		! Testing of the detailed balance
		!
		if (bal_test) then
		
			print *, "Detailed balance test"
			Ntest = 300
			mu    = 1
			nu    = 1
			
			w(1) = 10.0_dp*abs((Om(mu,nu))) 
			
			open(9,file="yy.dat")
			do i = 1, Nt
				yy(i) = &
				 exp(-ggA(i,mu) - (0.0,1.0)*llA(mu)*tt(i) &
				 	 -ggD(i,nu) - (0.0,1.0)*llD(nu)*tt(i)   )
				write(9,*) tt(i), real(yy(i)), aimag(yy(i)) 
			end do
			
			close(9)
		
			open(10,file="spectdens.dat")
			open(11,file="spectmod.dat")
			if (w(1) > 0) then
	    		wi = w(1)
	    		wc = w(1)/real(Ntest)
				do while(w(1) > -wi)  
					!wk(1) = w(1)
					call primitive(tt,real(yy),w,pyr)
					call primitive(tt,aimag(yy),w,pyi)
					integ = real(pyr(Nt,1)+(0.0_dp,1.0_dp)*pyi(Nt,1))
					write(10,'(2'//FMT_D//')') w(1), integ
					if (w(1) > 0.0_dp) then
						bal = exp(-((w(1)))*convFT/convET)
						write(11,'(2'//FMT_D//')') -w(1), integ*bal
					end if
					w(1) = w(1) - wc
				end	do
			end if
			close(10)
			close(11)
	
		end if
	
	end subroutine f_init_1
	
	!
	! Cleans Foerster module
	!
	subroutine clean_foerster()
		deallocate(foerster_rate,foerster_overlap)
	end subroutine clean_foerster

end module sci_foerster
