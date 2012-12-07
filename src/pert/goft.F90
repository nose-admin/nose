#include "util_allocation.h"

module goft

    use excitons

    implicit none

    ! Brownian
    character(len=32), parameter :: MT_0 = "BROWNIAN"                   ! implemented
    character(len=32), parameter :: MT_1 = "BROWNIAN_GENERAL"
    character(len=32), parameter :: MT_2 = "BROWNIAN_OVERDUMPED"
    character(len=32), parameter :: MT_3 = "BROWNIAN_UNDERDUMPED"

    ! Delta correlated
    character(len=32), parameter :: MT_4 = "DELTA"                      ! implemented

    ! Generic models
    character(len=32), parameter :: MT_10 = "OHMIC"
    character(len=32), parameter :: MT_11 = "JANG_CAO_SILBEY"
    character(len=32), parameter :: MT_12 = "DEBYE" ! same as overdumped BO?

    ! Special models
    character(len=32), parameter :: MT_51 = "RENGER_BCHL"

    type pag
        complex(dpc), dimension(:), pointer :: gg => NULL()
    end type

    type(pag), dimension(:), pointer :: all_goft => NULL()
    type(pag), dimension(:), pointer :: all_coft => NULL()
    type(pag), dimension(:), pointer :: all_hoft => NULL()

    !private :: pag, all_goft, all_coft
    private::brownian_general
    private::brownian
    private::brownian_underdamped
    private::deltacf
    private::dimensionless_CC
    public::ctanh

	logical :: g21_allocated, g2_allocated


contains

    !
    ! Initialization of the c(t), h(t) and g(t) functions
    !
    subroutine init_goft(err)
        integer, intent(out) :: err

        character(len=256) :: cbuff

        err = 0

		g21_allocated = .false.
		g2_allocated = .false.

        ! go through all gofts and create their time dependecies
        call resources_rewind_gofts()

        ! create a pointer to all time dependencies
        if (.not.associated(all_goft)) then
        	ALLOCATE(all_goft,(nr_gofts))
        	ALLOCATE(all_coft,(nr_gofts))
        	ALLOCATE(all_hoft,(nr_gofts))
        end if

        do

            write(cbuff,'(a,i3)') "Making g(t) - id : ",current_s_goft%id
            call print_log_message(trim(cbuff),6)
            call goft_create(current_s_goft%types, &
                             current_s_goft%params, &
                             current_s_goft%gt, &
                             current_s_goft%ct, &
                             current_s_goft%ht, &
                             current_s_goft%lambda)


            all_goft(current_s_goft%id)%gg => current_s_goft%gt
            all_hoft(current_s_goft%id)%gg => current_s_goft%ht
            all_coft(current_s_goft%id)%gg => current_s_goft%ct


            if (save_goft=="yes") then
                call goft_save_time_dependence(current_s_goft)
                call hoft_save_time_dependence(current_s_goft)
                call coft_save_time_dependence(current_s_goft)
            end if

            if (.not.resources_have_next_goft()) exit

            call resources_next_goft()

        end do

       ! set all_goft

    end subroutine init_goft




    !
    ! Calculates g(t) of the eigenstates
    !
    subroutine goft_prepare_excitonic_goft()
        integer :: i,j,k,a, b,c,d
        integer :: ij
        real(dp) :: cmpg, cmpgs, cmpgs2, ccount

        logical :: diag_only

        real(dp) :: rab

        real(dp), dimension(:,:), allocatable :: ST_2, STL_2, STR_2

		call print_log_message("Preparing excitonic goft ...", 8)

        call resources_rewind_blocks()

        do

            ! Reorganization energy of i-th eigenstate
            current_e_block%ll = 0.0_dp
            do i = 1,N1
                do j = 1,N1
                    !print *, i,j,N1,current_s_block%gindex(j),current_s_block%ll(i)
                     current_e_block%ll(i) =  current_e_block%ll(i) + &
                        (current_e_block%SS(j,i)**4)* &
                        current_s_block%ll(j) !current_s_block%gindex(j))
                end do
            end do


            cmpg = 0.0_dp
            cmpgs = 0.0_dp
            cmpgs2 = 0.0_dp
            ccount = 0

            diag_only = .true.

            if (diag_only) then

            do k = 1, Nte
                do i = 1,N1
                	rab = 0

                        j = i
                        ij = i !(i-1)*N1 + j
                        current_e_block%gg(ij,k) = 0.0_dp
                        current_e_block%cc(ij,k) = 0.0_dp
                        do a = 1, N1

                            !
                            ! here we save the excitonic g(t) function
                            !
                            current_e_block%gg(ij,k) = current_e_block%gg(ij,k) + &
                                (current_e_block%SS(a,i)**2)*(current_e_block%SS(a,j)**2) &
                                *all_goft(current_s_block%gindex(a))%gg(1+gt(1)*(k-1))
                            current_e_block%cc(ij,k) = current_e_block%cc(ij,k) + &
                                (current_e_block%SS(a,i)**2)*(current_e_block%SS(a,j)**2) &
                                *all_coft(current_s_block%gindex(a))%gg(1+gt(1)*(k-1))
                            if (k == 1) then
                                cmpg = cmpg + (current_e_block%SS(a,i)**2)*(current_e_block%SS(a,j)**2)
                            end if

                        end do

                        if (k == 1) then
                        cmpgs = cmpgs + cmpg
                        if (cmpg > 1.0/real(N1)) then
                                cmpgs2 = cmpgs2 + cmpg
                                ccount = ccount + 1
                        end if

                        cmpg = 0.0_dp
                        if (j == N1) then
                            cmpgs = 0.0_dp
                            cmpgs2 = 0.0_dp
                            ccount = 0
                        end if
                        end if

                end do
            end do

            else

            do k = 1, Nte
                do i = 1,N1
                    do j = 1,N1
                        ij = (i-1)*N1 + j
                        current_e_block%gg(ij,k) = 0.0_dp
                        do a = 1, N1

                            !
                            ! here we save the excitonic g(t) function
                            !
                            current_e_block%gg(ij,k) = current_e_block%gg(ij,k) + &
                                (current_e_block%SS(a,i)**2)*(current_e_block%SS(a,j)**2) &
                                *all_goft(current_s_block%gindex(a))%gg(1+gt(1)*(k-1))
                            if (k == 1) then
                                cmpg = cmpg + (current_e_block%SS(a,i)**2)*(current_e_block%SS(a,j)**2)
                            end if
                        end do
                        if (k == 1) then
                        cmpgs = cmpgs + cmpg
                        if (cmpg > 1.0/real(N1)) then
                                cmpgs2 = cmpgs2 + cmpg
                                ccount = ccount + 1
                        end if
                        write(55,*) i,j,cmpg
                        cmpg = 0.0_dp
                        if (j == N1) then
                            write(55,*) "Suma   ", cmpgs
                            write(55,*) "Suma 2 ", cmpgs2, ccount
                            cmpgs = 0.0_dp
                            cmpgs2 = 0.0_dp
                            ccount = 0
                        end if
                        end if
                    end do
                end do
            end do

            end if


            ! correlation functions of two-exciton states

			if (use_twoexcitons) then

				! association status and size of the unused pointer seem to be different onder different compilations
				! with or without mpi
!				if ((.not.associated(current_e_block%gg_2)).or.(size(current_e_block%gg_2)<=1)) then
				if (.not.g2_allocated) then

					ALLOCATE(current_e_block%gg_2, (N1*(N1-1)/2,N1*(N1-1)/2,Nte))
					g2_allocated = .true.

				end if

!				if ((.not.associated(current_e_block%gg_21)).or.(size(current_e_block%gg_21)<=1)) then
				if (.not.g21_allocated) then

					ALLOCATE(current_e_block%gg_21,(N1*(N1-1)/2,N1,Nte))
					g21_allocated = .true.

				end if

				ALLOCATE(ST_2, (N1,N1))
				ALLOCATE(STL_2, (N1*(N1-1)/2,N1*(N1-1)/2) )
				ALLOCATE(STR_2, (N1*(N1-1)/2,N1*(N1-1)/2) )

				call init_kroneker(N1)

				ST_2  = abs(current_e_block%SS)**2
				STL_2 = abs(current_e_block%SS_2)**2
				STR_2 = transpose(STL_2)

            	do k = 1, Nte


            	    do j = 1, N1*(N1-1)/2


                		do i = 1,N1

                		a = current_e_block%ione1(j)
                		b = current_e_block%ione2(j)

                    	current_e_block%gg_21(j,i,k) = current_e_block%gg(a,k)*kroneker(a,i)  &
                    								 +  current_e_block%gg(b,k)*kroneker(b,i)


                    	end do

                    	do i = 1,N1*(N1-1)/2

                    	c = current_e_block%ione1(i)
                		d = current_e_block%ione2(i)


                    	current_e_block%gg_2(j,i,k) = current_e_block%gg(a,k)*(kroneker(a,c)+kroneker(a,d)) &
                    								+ current_e_block%gg(b,k)*(kroneker(b,c)+kroneker(b,d))


                    	end do

                    end do

					! Transform into the excitonic basis

					current_e_block%gg_21(:,:,k) = matmul(STL_2,matmul(current_e_block%gg_21(:,:,k),ST_2))


					current_e_block%gg_2(:,:,k) = matmul(STL_2,matmul(current_e_block%gg_2(:,:,k),STR_2))

                 end do

                 call delete_kroneker()


			end if

            if (.not.resources_have_next_block()) exit

            call resources_next_block()


        end do

		DEALLOCATE(ST_2)
	    DEALLOCATE(STL_2)
		DEALLOCATE(STR_2)


		call print_log_message("... excitonic goft done.", 8)


    end subroutine goft_prepare_excitonic_goft

    !
    ! Calculate the time dependence of a given g(t) function
    !
    ! tp ...... types of the g(t) modes
    ! params .. parameters of the modes
    ! gt ...... time dependence
    !
    subroutine goft_create(tp,params,ggt,cct,hht,lambda)
        character(len=*), dimension(:), intent(in) :: tp
        real(dp), dimension(:,:), intent(in)       :: params
        real(dp), intent(out)                      :: lambda
        complex(dpc), dimension(:), pointer        :: ggt
        complex(dpc), dimension(:), pointer        :: cct
        complex(dpc), dimension(:), pointer        :: hht
        !local
        integer :: i, j
        real(dp) :: ll

        lambda = 0.0_dp
        ll = 0.0_dp

        if (.not.associated(ggt)) then
            ALLOCATE(ggt,(grid_Nt))
            ALLOCATE(cct,(grid_Nt))
            ALLOCATE(hht,(grid_Nt))
            ggt = 0.0_dp
            cct = 0.0_dp
            hht = 0.0_dp
        else
            print *, size(ggt,1), grid_Nt
            if (size(ggt,1) /= grid_Nt) stop "Wrong size of the g(t) vector"
        end if


        do i = 1,(size(tp,1))

            if (trim(tp(i)) == "BROWNIAN") then

                call brownian(params(1:2,i),ggt,cct,hht,ll,ADD="yes")
                lambda = lambda + ll

            else if (trim(tp(i)) == "DELTA") then

                !print *, "GAMMA = ", params(1,i)
                call deltacf(params(1,i),ggt,cct,hht,ADD="yes")

            else if (trim(tp(i)) == "BROWNIAN_GENERAL") then

                call brownian_general(params(1:3,i),ggt,cct,hht,ll,ADD="yes")
                lambda = lambda + ll

            else if (trim(tp(i)) == "BROWNIAN_UNDERDAMPED") then

                call brownian_underdamped(params(1:2,i),ggt,cct,hht,ll,ADD="yes")
                lambda = lambda + ll

			else if (trim(tp(i)) == "GAUSSIAN") then

				call gaussiancf(params(1:2,i),ggt,cct,hht,ll,ADD="yes")
                lambda = lambda + ll

            end if

        end do



    end subroutine goft_create

    !
    ! Saves the time dependence of a given line shape function
    !
    subroutine goft_save_time_dependence(goft)
        type(site_goft), intent(in) :: goft

        ! local
        integer :: i
        character(len=10) :: cbuff
        complex(dpc), dimension(:), pointer :: gg => NULL()

        gg => goft%gt

        write(cbuff,'(i1)') goft%index
        open(unit=10,file=trim(file_join("log","goft"//trim(cbuff)//".dat")))

        do i = 1, grid_Nt
            write(10,*) grid_t(i), real(gg(i)), aimag(gg(i))
        end do

        close(10)

    end subroutine goft_save_time_dependence

    !
    ! Saves the time dependence of a given line shape function
    !
    subroutine hoft_save_time_dependence(goft)
        type(site_goft), intent(in) :: goft

        ! local
        integer :: i
        character(len=10) :: cbuff
        complex(dpc), dimension(:), pointer :: gg => NULL()

        gg => goft%ht

        write(cbuff,'(i1)') goft%index
        open(unit=10,file=trim(file_join("log","hoft"//trim(cbuff)//".dat")))

        do i = 1, grid_Nt
            write(10,*) grid_t(i), real(gg(i)), aimag(gg(i))
        end do

        close(10)

    end subroutine hoft_save_time_dependence

    !
    ! Saves the time dependence of a given line shape function
    !
    subroutine coft_save_time_dependence(goft)
        type(site_goft), intent(in) :: goft

        ! local
        integer :: i
        character(len=10) :: cbuff
        complex(dpc), dimension(:), pointer :: gg => NULL()

        gg => goft%ct

        write(cbuff,'(i1)') goft%index
        open(unit=10,file=trim(file_join("log","coft"//trim(cbuff)//".dat")))

        do i = 1, grid_Nt
            write(10,*) grid_t(i), real(gg(i)), aimag(gg(i))
        end do

        close(10)

    end subroutine coft_save_time_dependence


    !
    ! Brownian mod of g(t)
    !
    ! This function is obsolete and should be removed in future. It calculates
    ! the overdumped HO CF with coth instead of cot. (Differs from second order
    ! in beta-expansion.
    !
    !
    subroutine brownian_old(params,ggt,cct,hht,lambda,ADD)
        real(dp), dimension(:), intent(in)       :: params
        complex(dpc), dimension(:), pointer      :: ggt
        complex(dpc), dimension(:), pointer      :: cct,hht
        real(dp), intent(out)                    :: lambda
        character(len=*), intent(in), optional   :: ADD
        ! local
        integer  :: i,j,k  ! indices
        real(dp) :: t    ! time
        integer  :: Ntt  ! number of time steps
        integer  :: Nmatsu ! number of Matsubara terms

        real(dp) :: C2r,C2i,C1r,C1i,Cr,Ci,C2rq,C2iq,C2r_loc,C2i_loc
        real(dp) :: ssum, ssum1, ssum2, ssum3, rpom, diff
        real(dp) :: BH, mfac

        real(dp), dimension(1) :: Kfac, Sfac,Sfac1,Sfac2,RR,ll,Cfac,eps,inhomo
        real(dp), dimension(0:100,1) :: cc,c1,c2

        !
        ! Set evaluation parameters
        !
        ll(1) = params(1) !Energy_cm_to_internal*params(1)
        lambda = params(1)
        RR(1) = 1.0_dp/params(2)
        eps = 1.0_dp
        inhomo = 0.0_dp

        ! dt ... elementary grid step

        Ntt = grid_Nt

        BH = (0.6582120_dp/8.617385d-5)/temp  ! hbar/kT
        !print *, "Temparature: ", temp

        ! estimation fo the Matsubara contribution for a given temperature

        ! matsubara factor
        mfac = 2.0_dp*PI_D/BH
        cc = 0.0_dp

        !do i = 1, Nc
        i = 1

!            if ((mtype(i) == MT_2).or.(mtype(i) == MT_0)) then

                ! another factor that appeares in the sum
                Kfac(i) = (RR(i)*BH/(2.0_dp*PI_D))**2
                ! sum prefactor
                Sfac(i) = 2.0_dp*RR(i)*ll(i)/PI_D
                Sfac1(i) = ll(i)*RR(i)*BH/(PI_D**2)
                Sfac2(i) = 4.0_dp*ll(i)*RR(i)*(BH**2)/((2.0_dp*PI_D)**3)

                ! coth prefactor
                Cfac(i) = ll(i)*RR(i)*(1.0_dp/tanh(BH*RR(i)/2.0_dp))

                Nmatsu = 0
                ssum = 0.0_dp
                do

                    Nmatsu = Nmatsu + 1
                    if (Nmatsu > 100) exit

                    cc(Nmatsu,i) = real(Nmatsu,dp)/(real(Nmatsu**2,dp)-Kfac(i))
                    !print *, Nmatsu, cc(Nmatsu,i)
                    c1(Nmatsu,i) = 1.0_dp/(real(Nmatsu**2,dp)-Kfac(i))
                    c2(Nmatsu,i) = 1.0_dp/(real(Nmatsu,dp)*(real(Nmatsu**2,dp)-Kfac(i)))
                    ssum = ssum + cc(Nmatsu,i)
                    diff = abs((cc(Nmatsu,i)-cc(Nmatsu-1,i))/ssum)*100
                    !print *, Nmatsu, diff
                    !if (diff < 1.0_dp) then
                    if (Nmatsu == 20) then
                        !print *, "convergency: ", diff
                        diff = abs(Sfac(i)*ssum/Cfac(i))
                        !print *, "Nmatsu = ",Nmatsu, " ratio: ", diff
                        exit
                    end if

                end do


!            end if


        !end do
!		open(unit=11,file='/home/olsij4am/prace/NOSE-DEBUG.dat')
!		open(unit=12,file='/home/olsij4am/prace/NOSE-DEBUG2.dat')

        do i = 1, Ntt

            t = grid_t(i)

            C2r = 0.0_dp
            C2i = 0.0_dp
            C1r = 0.0_dp
            C1i = 0.0_dp
            Cr  = 0.0_dp
            Ci  = 0.0_dp

            C2rq = 0.0_dp
            C2iq = 0.0_dp

            !tssum = 0.0_dp

            ! sum over modes
            !do j = 1, Nc
            j = 1

            ! Brownian oscillator g(t),g'(t),g''(t)
!            if ((mtype(j) == MT_2).or.(mtype(j)== MT_0)) then

                ! Matsubara sum
                ssum = 0.0_dp
                ssum1 = 0.0_dp
                ssum2 = 0.0_dp
                ssum3 = 0.0_dp
                do k = 1, Nmatsu
                    ssum  = ssum  + cc(k,j)*(exp(-real(k,dp)*mfac*t))
                    ssum1 = ssum1 + c1(k,j)*(exp(-real(k,dp)*mfac*t)-1.0_dp)
                    ssum2 = ssum2 + c2(k,j)*(exp(-real(k,dp)*mfac*t)-1.0_dp)
                    ssum3 = ssum3 + c1(k,j)*t
                end do
                !tssum  = tssum  + Sfac(j)*ssum
                !tssum1 = tssum1 + Sfac1(j)*ssum1
                !tssum2 = tssum2 + Sfac2(j)*ssum2 - Sfac1(j)*ssum3

                rpom = 1.0_dp/tanh(BH*RR(j)/2.0_dp)

                C2r_loc = (ll(j)*RR(j)*rpom*exp(-RR(j)*t) + Sfac(j)*ssum)
                C2i_loc = - (ll(j)*RR(j)*exp(-RR(j)*t))

                ! g''(t)
                C2r = C2r + C2r_loc/(eps(j)**2)
                C2i = C2i + C2i_loc/(eps(j)**2)

                C2rq = C2rq + ((((1-eps(j)**2)/eps(j)**2)**2)/8)*(C2r_loc/ll(j))**2
                C2iq = C2iq + ((((1-eps(j)**2)/eps(j)**2)**2)/8)*(C2i_loc/ll(j))**2

                ! g'(t)
                C1r = C1r - ll(j)*rpom*(exp(-RR(j)*t) - 1.0_dp) - Sfac1(j)*ssum1
                C1i = C1i + ll(j)*(exp(-RR(j)*t) - 1.0_dp)

                ! g(t)
                Cr  = Cr  + (ll(j)/RR(j))*rpom*(exp(-RR(j)*t)+RR(j)*t-1.0_dp)  &
                            + Sfac2(j)*ssum2 + Sfac1(j)*ssum3 + &
                            ((inhomo(j)**2)/2.0_dp)*(t**2)
                Ci  = Ci  - (ll(j)/RR(j))*(exp(-RR(j)*t)+RR(j)*t-1.0_dp)

 !           end if


        !end do

!        write(11,*) i*dt,real(C2r + (0.0_dp,1.0_dp)*C2i),real(C1r + (0.0_dp,1.0_dp)*C1i),real(Cr + (0.0_dp,1.0_dp)*Ci)
!        write(12,*) i*dt,aimag(C2r + (0.0_dp,1.0_dp)*C2i),aimag(C1r + (0.0_dp,1.0_dp)*C1i),aimag(Cr + (0.0_dp,1.0_dp)*Ci)


        if (present(ADD)) then

            ggt(i) = ggt(i) + Cr + (0.0_dp,1.0_dp)*Ci
            hht(i) = hht(i) + C1r + (0.0_dp,1.0_dp)*C1i
            cct(i) = cct(i) + C2r + (0.0_dp,1.0_dp)*C2i

        else

            ggt(i) = Cr + (0.0_dp,1.0_dp)*Ci
            hht(i) = C1r + (0.0_dp,1.0_dp)*C1i
            cct(i) = C2r + (0.0_dp,1.0_dp)*C2i

        end if


    end do

!    close(11)
!    close(12)
!    write(*,*) 'debug functions written, ending program'
    !stop

    end subroutine brownian_old


    !
    ! Complex tanh
    !
	function ctanh(x) result(nav)
		complex(dpc), intent(in)	:: x
		complex(dpc)				:: nav

		nav = (exp(x)-exp(-x))/(exp(x)+exp(-x))
	end function ctanh


    !
    ! Dimensionless brownian overdumped CC(T,x),
    !    x = beta hbar Lambda/2
    !    T = 2 t / beta hbar
    !    C(Lambda, ll, t) = Lambda ll CC(T,x)
    !
    !    CC(T,x) = (cot(x)-i) e^(-Tx) + sum_n=1^\infty e^(-n pi T) (1/(x + n pi) - 1/(x - n pi))
    !
	function dimensionless_CC(T,x) result(CC)
		real(dp), intent(in)	    :: x, T
		complex(dpc)				:: CC

		real(dp)					:: y, diff
		integer(i4b) 				:: nearest_pole,i
		character(len=256)			:: buff

		! main idea of this function is to evaluate smooth CC(T,x) by
		! proper handling of summation of the Matsubara sum by switching
		! between Laurent expansion around singularities of cot(x),
		! which are regularized by proper terms of the sum, and
		! evaluating cot(x) exactly outside the singularities

		! we shift x to the nearest singularity
		nearest_pole = INT((x + PI_D/2) / (PI_D))
		y = x - nearest_pole*PI_D

		! if assures that the message is written only once - at time = 0
		if(T < 1e-5) then
			write(buff,'(a,f6.3,a)') 'brownian cf taken in ',x/PI_D,', where 1 = nearest singularity'

			if(x/PI_D < 0.8) then
				call print_log_message(trim(buff),7)
			else
				call print_warning_message(trim(buff),5)
			end if
		end if

		! if we are +-1 around pole, Laurent expansion works well, otherwise
		! we are far enough and we don't need it
		if(abs(y) < 0.1 .and. nearest_pole > 0) then

			if(T < 1e-5) then
				call print_warning_message('interpolating cf around pole',5)
			end if

			! series of e^(-T*y) cot y - 1/y
			CC = 0
			CC = CC + (y**3) * (15*(T**4) - 60*(T**2) - 8.0_dp)/360.0_dp
			CC = CC + (y**2) * (T/3 - (T**3)/6)
			CC = CC + (y**1) * (-1.0_dp/3.0_dp + T**2/2.0_dp)
			CC = CC - T

			! second partial fraction (+1/(y + 2 pi n)) is added exactly
			CC = CC + 1.0_dp/(y + 2.0_dp*PI_D*nearest_pole)

			! imaginary part is added exactly
			CC = CC - exp(-T*y)*cmplx(0,1,dpc)

			! whole previous CC is multiplied by exp(-pi n T)
			CC = CC * exp(- PI_D*nearest_pole*T)

		else

			CC =      exp(-T*x)*(1.0_dp/tan(x) - cmplx(0,1,dpc))

			if(nearest_pole > 0) then
				CC = CC + exp(- PI_D*nearest_pole*T)*						&
							(1.0_dp/(x+nearest_pole*PI_D)-1.0_dp/(x-nearest_pole*PI_D))
			end if
		end if

		! the rest of Matsubara sum is added, except the nearest pole
		do i = 1, 2000
			if(i == nearest_pole) then
				cycle
			end if

			diff = exp(- PI_D*i*T)*(1.0_dp/(x+i*PI_D)-1.0_dp/(x-i*PI_D))
			CC = CC + diff

			if(diff/abs(CC) < 0.001_dp) then
				exit
			end if
		end do

		if(T < 1e-5) then
			write(buff,'(a,i4)') 'Matsubara sum to ',i-1
			call print_log_message(trim(buff),7)
		end if


	end function dimensionless_CC

	subroutine brownian(params,ggt,cct,hht,lambda,ADD)
        real(dp), dimension(:), intent(in)		:: params
        complex(dpc), dimension(:), pointer		:: ggt
        complex(dpc), dimension(:), pointer		:: cct,hht
        real(dp), intent(out)						:: lambda
        character(len=*), intent(in), optional	:: ADD

        complex(dpc), dimension(size(ggt))	:: cct_tmp,hht_tmp,ggt_tmp
        real(dp)								:: BH, LLambda,t
        integer(i4b)							:: Ntt, i

        !
        ! Set evaluation parameters
        !
        lambda 	= params(1)
        LLambda 	= 1.0_dp/params(2)

        ! dt ... elementary grid step

        Ntt = size(ggt)

        BH = (0.6582120_dp/8.617385d-5)/temp  ! hbar/kT

        do i=1, Ntt
        	t = (i-1)*dt

	        cct_tmp(i) = lambda*LLambda*dimensionless_CC(2.0*t/BH, BH*LLambda/2.0)
	        if(i > 1) then
    	    	hht_tmp(i) = hht_tmp(i-1) + dt*cct_tmp(i)
        		ggt_tmp(i) = ggt_tmp(i-1) + dt*hht_tmp(i)
        	else
    	    	hht_tmp(i) = dt*cct_tmp(i)
        		ggt_tmp(i) = dt*hht_tmp(i)
        	end if
        end do

!		open(unit=11,file='/home/olsij4am/prace/NOSE-debug.dat')
!		open(unit=12,file='/home/olsij4am/prace/NOSE-debug2.dat')

    	! write to global functions
        if (.not. present(ADD)) then
   	    	ggt = 0.0_dp
       		cct = 0.0_dp
       		hht = 0.0_dp
        end if

	   	do i=1,Ntt
 !       	write(11,*) i*dt,real(ggt_tmp(i)),real(hht_tmp(i)),real(cct_tmp(i))
  !      	write(12,*) i*dt,aimag(ggt_tmp(i)),aimag(hht_tmp(i)),aimag(cct_tmp(i))

        	cct(i) = cct_tmp(i) + cct(i)
    	    hht(i) = hht_tmp(i) + hht(i)
        	ggt(i) = ggt_tmp(i) + ggt(i)
    	end do

!	   	close(11)
!   	close(12)
!    	write(*,*) 'debug functions written'
!    	stop


	end subroutine brownian


	subroutine brownian_underdamped(params,ggt,cct,hht,lambda,ADD)
        real(dp), dimension(:), intent(in)		:: params
        complex(dpc), dimension(:), pointer		:: ggt
        complex(dpc), dimension(:), pointer		:: cct,hht
        real(dp), intent(out)						:: lambda
        character(len=*), intent(in), optional	:: ADD

        complex(dpc), dimension(size(ggt))	:: cct_tmp,hht_tmp,ggt_tmp
        real(dp)								:: BH, t, omega
        complex(dpc)							:: hht_tmp0, ggt_tmp0
        integer(i4b)							:: Ntt, i

        !
        ! Set evaluation parameters
        !
        lambda 	= params(1)
        omega = params(2)

        ! dt ... elementary grid step

        Ntt = size(ggt)

        BH = (0.6582120_dp/8.617385d-5)/temp  ! hbar/kT

        do i=1, Ntt
        	t = (i-1)*dt

!			! rozmer casu z LLambda, ktera je jinak predfaktorem -- sedi?
!	        cct_tmp(i) = lambda*omega*(cos(t*omega)/tanh(omega/2*BH)-CMPLX(0.0_dp,1.0_dp)*sin(t*omega))
!	        if(i > 1) then
!    	    	hht_tmp(i) = hht_tmp(i-1) + dt*cct_tmp(i)
!        		ggt_tmp(i) = ggt_tmp(i-1) + dt*hht_tmp(i)
!        	else
!    	    	hht_tmp(i) = dt*cct_tmp(i)
!        		ggt_tmp(i) = dt*hht_tmp(i)
!        	end if

			! ANALYTICAL
        	hht_tmp(i) = lambda*(sin(t*omega)/tanh(omega/2*BH)-CMPLX(0.0_dp,1.0_dp)*(-cos(t*omega)))

        	if(i == 1) then
        		hht_tmp0 = hht_tmp(1)
        	end if

        	hht_tmp(i) = hht_tmp(i) - hht_tmp0

        	ggt_tmp(i) = lambda/omega*(-cos(t*omega)/tanh(omega/2*BH)-CMPLX(0.0_dp,1.0_dp)*(-sin(t*omega)))

        	if(i == 1) then
        		ggt_tmp0 = ggt_tmp(1)
        	end if

        	ggt_tmp(i) = ggt_tmp(i) - ggt_tmp0
        end do

!		open(unit=11,file='/home/olsij4am/prace/NOSE-debugU.dat')
!		open(unit=12,file='/home/olsij4am/prace/NOSE-debugU2.dat')

    	! write to global functions
        if (.not. present(ADD)) then
   	    	ggt = 0.0_dp
       		cct = 0.0_dp
       		hht = 0.0_dp
        end if

	   	do i=1,Ntt
!			write(11,*) i*dt,real(ggt_tmp(i)),real(hht_tmp(i)),real(cct_tmp(i))
!			write(12,*) i*dt,aimag(ggt_tmp(i)),aimag(hht_tmp(i)),aimag(cct_tmp(i))

        	cct(i) = cct_tmp(i) + cct(i)
    	    hht(i) = hht_tmp(i) + hht(i)
        	ggt(i) = ggt_tmp(i) + ggt(i)
    	end do

!	   	close(11)
!	   	close(12)
!	   	write(*,*) 'debug functions written'
!		write(*,*) 'lambda, omega', lambda, omega
!    	stop


	end subroutine brownian_underdamped


    !
    ! General brownian mod of g(t)
    !
    subroutine brownian_general(params,ggt,cct,hht,lambda,ADD)
		real(dp), dimension(:), intent(in)       :: params
        complex(dpc), dimension(:), pointer      :: ggt
        complex(dpc), dimension(:), pointer      :: cct,hht

        complex(dpc), dimension(size(ggt))      :: cct_tmp,hht_tmp,ggt_tmp
        real(dp), intent(out)                    :: lambda
        character(len=*), intent(in), optional   :: ADD
        ! local
        integer  :: i,j,k  ! indices
        real(dp) :: t    ! time
        integer  :: Ntt  ! number of time steps
        integer  :: Nmatsu ! number of Matsubara terms

        real(dp) :: omega,gamma,ll
        real(dp) :: C2r,C2i
        real(dp) :: ssum
        real(dp) :: BH, mfac, xi_square_divided_by_mass
        real(dp), dimension(1) :: www

        !real(dp), dimension(2) :: debug_params

        complex(dp) :: Phi,Phi_prime,zeta
        real(dp), dimension(0:100) :: cc

		if(.not. (size(params) == 3)) then
        	call print_error_message(-1,"wrong number of params in brownian_general() (goft.F90)")
        end if

        !! debug - MUST be removed when calculating
        !debug_params(1) = params(1)
        !debug_params(2) = params(3)
        !call brownian_underdamped(debug_params(1:2),ggt,cct,hht,lambda,ADD)

        !
        ! Set evaluation parameters
        !
        ll = params(1) !Energy_cm_to_internal*params(1)
        lambda = params(1)
        gamma = 1.0_dp/params(2)
        omega = params(3)


        zeta 		= sqrt(cmplx(omega**2 - (gamma**2)/4.0_dp,0.0_dp,dpc))
        Phi  		= gamma/2.0_dp + cmplx(0,1,dpc)*zeta
        Phi_prime	= gamma/2.0_dp - cmplx(0,1,dpc)*zeta

        xi_square_divided_by_mass = lambda*(gamma**2+4*zeta**2)/2

        if(real(zeta) < 0) then
        	call print_error_message(-1,"sqrt with negative real part in goft!")
        end if

        Ntt = size(ggt)

        BH = (0.6582120_dp/8.617385d-5)/temp  ! hbar/kT

        ! matsubara factor
        mfac = 2.0_dp*PI_D/BH
        cc = 0.0_dp

        Nmatsu = 0
        do
        	Nmatsu = Nmatsu + 1

	        cc(Nmatsu) = 2*gamma/BH *   &
	        				Nmatsu*mfac/( (( (omega**2) + ((Nmatsu*mfac)**2) )**2) - (gamma**2) * ((mfac*Nmatsu)**2) )

            if (Nmatsu >= 100) then
            	exit
            end if
        end do

        do i = 1, Ntt

            t = grid_t(i)

            C2r = 0.0_dp
            C2i = 0.0_dp

            j = 1

                ! Matsubara sum
                ssum = 0.0_dp
                do k = 1, Nmatsu
                    ssum  = ssum  + cc(k)*(exp(-real(k,dp)*mfac*t))
                end do

				! g''(t)
                C2r =  	1.0_dp/ctanh(cmplx(0,1,dpc)*BH*Phi_prime/2.0_dp)*exp(-Phi_prime*t) - &
                		1.0_dp/ctanh(cmplx(0,1,dpc)*BH*Phi/2.0_dp)*exp(-Phi*t)
				C2r = 	C2r/4.0_dp/zeta + ssum
                C2i = - 1.0_dp/2.0_dp/zeta * exp(-gamma*t/2.0_dp)*sin(zeta*t)

                C2r = C2r * xi_square_divided_by_mass
                C2i = C2i * xi_square_divided_by_mass

			! we integrate h(t), c(t)
			!! better integration should be applied

	        cct_tmp(i) = C2r + (0.0_dp,1.0_dp)*C2i

    	end do

    	www = 0.0_dp
    	call primitive_cmplx(grid_t, cct_tmp, www, hht_tmp)
    	call primitive_cmplx(grid_t, hht_tmp, www, ggt_tmp)

    	if(2*PI/omega/dt < 50) then
    		call print_warning_message("Step too small, correlation function may be integrated wrongly", -1)
    	end if

!!    	! normalization to reorganization energy
!!		do i=1,Ntt
!!			cct_tmp(i) = -cct_tmp(i)/aimag(hht_tmp(Ntt))*ll
!!			ggt_tmp(i) = -ggt_tmp(i)/aimag(hht_tmp(Ntt))*ll
!!			hht_tmp(i) = -hht_tmp(i)/aimag(hht_tmp(Ntt))*ll
!!		end do

!		open(unit=11,file='/home/olsij4am/prace/nose-debug.dat')
!		open(unit=12,file='/home/olsij4am/prace/nose-debug2.dat')

    	! write to global functions
        if (.not. present(ADD)) then
   	    	ggt = 0.0_dp
       		cct = 0.0_dp
       		hht = 0.0_dp
        end if

    	do i=1,Ntt
!	       	write(11,*) i*dt,real(ggt_tmp(i)),real(hht_tmp(i)),real(cct_tmp(i))
!	       	write(12,*) i*dt,aimag(ggt_tmp(i)),aimag(hht_tmp(i)),aimag(cct_tmp(i))

        	cct(i) = cct_tmp(i) + cct(i)
    	    hht(i) = hht_tmp(i) + hht(i)
        	ggt(i) = ggt_tmp(i) + ggt(i)
    	end do

!		close(11)
!		close(12)
!		write(*,*) 'debug functions written'
!		stop

    end subroutine brownian_general


    !
    ! Delta correlated mod of g(t)
    !
    subroutine deltacf(params,ggt,cct,hht,ADD)
        real(dp), intent(in)       :: params
        complex(dpc), dimension(:), pointer      :: ggt
        complex(dpc), dimension(:), pointer      :: cct,hht
        character(len=*), intent(in), optional   :: ADD

        integer :: i, Ntt
        real(dp) :: t    ! time

        Ntt = grid_Nt

        do i = 1, Ntt

            t = grid_t(i)

            if (present(ADD)) then

                ggt(i) = ggt(i) + params*t
                hht(i) = hht(i) + params
                !cct(i) = cct(i) + 0.0_dp

            else

                ggt(i) = params*t
                hht(i) = params
                cct(i) = 0.0_dp

            end if

            !print *, t, ggt(i)

        end do

    end subroutine deltacf


	subroutine gaussiancf(params,ggt,cct,hht,lambda,ADD)
        real(dp), dimension(:), intent(in)		:: params
        complex(dpc), dimension(:), pointer		:: ggt
        complex(dpc), dimension(:), pointer		:: cct,hht
        real(dp), intent(out)						:: lambda
        character(len=*), intent(in), optional	:: ADD

        complex(dpc), dimension(size(ggt))	:: cct_tmp,hht_tmp,ggt_tmp
        real(dp)								:: BH, t, tc
        complex(dpc)							:: hht_tmp0, ggt_tmp0
        integer(i4b)							:: Ntt, i

        !
        ! Set evaluation parameters
        !
        lambda 	= params(1)
        tc = params(2)

        ! dt ... elementary grid step

        Ntt = size(ggt)

        BH = (0.6582120_dp/8.617385d-5)/temp  ! hbar/kT

        do i=1, Ntt
        	t = (i-1)*dt

			! ANALYTICAL
			cct_tmp(i) = 2.0_dp*lambda*((1.0_dp/BH)-(0.0d0,1.0d0)*t/(tc**2))*exp(-(t/tc)**2)

        	hht_tmp(i) = 2.0_dp*lambda*((sqrt(PI_D)/(2.0_dp*BH))*tc*erf(t/tc)-(0.0d0,1.0d0)*(1.0_dp-exp(-(t/tc)**2)))

        	if(i == 1) then
        		hht_tmp0 = hht_tmp(1)
        	end if

        	hht_tmp(i) = hht_tmp(i) - hht_tmp0

        	ggt_tmp(i) = lambda*((tc*sqrt(PI_D)/BH)*(t*erf(t/tc)+tc*exp(-(t/tc)**2)/sqrt(PI_D))- &
        	            (0.0d0,1.0d0)*(t-sqrt(PI_D)*tc*erf(t/tc)/2.0d0))

        	if(i == 1) then
        		ggt_tmp0 = ggt_tmp(1)
        	end if

        	ggt_tmp(i) = ggt_tmp(i) - ggt_tmp0
        end do

!		open(unit=11,file='/home/olsij4am/prace/NOSE-debugU.dat')
!		open(unit=12,file='/home/olsij4am/prace/NOSE-debugU2.dat')

    	! write to global functions
        if (.not. present(ADD)) then
   	    	ggt = 0.0_dp
       		cct = 0.0_dp
       		hht = 0.0_dp
        end if

	   	do i=1,Ntt
!			write(11,*) i*dt,real(ggt_tmp(i)),real(hht_tmp(i)),real(cct_tmp(i))
!			write(12,*) i*dt,aimag(ggt_tmp(i)),aimag(hht_tmp(i)),aimag(cct_tmp(i))

        	cct(i) = cct_tmp(i) + cct(i)
    	    hht(i) = hht_tmp(i) + hht(i)
        	ggt(i) = ggt_tmp(i) + ggt(i)
    	end do


!	   	close(11)
!	   	close(12)
!	   	write(*,*) 'debug functions written'
!		write(*,*) 'lambda, omega', lambda, tc

!    	stop


	end subroutine gaussiancf



	!
	! Overdamped Brownian oscillator correlation in frequency domain
	!
	function comega_overdamped_brownian(om,lam,tau,kbT) result (rr)
		real(dp), intent(in) :: om, lam, tau, kbT
	    real(dp) :: rr, Ll

	    real(dp) :: abom, cc

		Ll = 1/tau

	    if (om == 0.0_dp) then

	    	rr = 4.0_dp*PI_D*kbT*lam/Ll
	    	return

	    end if

	    abom = abs(om)
	    cc = 4.0_dp*PI_D*lam*abom*Ll/((1.0_dp-exp(-abom/kbT))*(abom**2 + Ll**2))

	    if (om < 0) then
	    	rr = cc*exp(-abom/kbT)
	    else
	    	rr = cc
	    end if

	end function comega_overdamped_brownian


	!
	! Overdamped Brownian oscillator correlation in frequency domain
	!
	function comega_general_brownian(om,lam,tau,omj,kbT) result (rr)
		real(dp), intent(in) :: om, lam, tau, omj, kbT
	    real(dp) :: rr, Ll

	    real(dp) :: abom, cc

		Ll = 1/tau

	    if (om == 0.0_dp) then

	    	rr = 4.0_dp*PI_D*kbT*lam/Ll
	    	return

	    end if

	    abom = abs(om)
	    cc = 4.0_dp*PI_D*lam*abom*Ll/((1.0_dp-exp(-abom/kbT))*( ((abom**2 - omj**2)/abom)**2 + Ll**2))

	    if (om < 0) then
	    	rr = cc*exp(-abom/kbT)
	    else
	    	rr = cc
	    end if

	end function comega_general_brownian




end module goft
