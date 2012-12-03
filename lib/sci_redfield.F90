!
! Redfield relaxation rates
!
module sci_redfield

   use numer_matrix
   use numer_interp
   use std_io

   implicit none

   real(dp), dimension(:,:), allocatable :: redfield_rates_1ex
!   real(dp), dimension(:,:,:,:), allocatable :: redfield_rate_2ex

   real(dp), dimension(:,:), allocatable :: redfield_deph_1ex
   real(dp), dimension(:),   allocatable :: redfield_deph_01ex

   interface init_redfield
   		module procedure init_redfield_1, init_redfield_all
   end interface init_redfield

   private :: init_redfield_1, init_redfield_all

contains


    subroutine init_redfield_1(en)
        real(dp), intent(in), dimension(:) :: en        ! energies

		integer :: N

        N = size(en,1)   ! number of energis

        !
        ! allocate only first time you are here
        !
        if (allocated(redfield_rates_1ex)) then
            if (size(redfield_rates_1ex,1)/=N) then
                call clean_redfield()
                call allocate_redfield(N)
            end if
        else
            call allocate_redfield(N)
        end if

 	end subroutine init_redfield_1


    !
    ! Initialization of the Redfield tensor
    !
    !
    subroutine init_redfield_all(en,ll,SS,cc,dt,Temp)
        real(dp), intent(in), dimension(:) 	:: en        ! energies
        real(dp), intent(in), dimension(:) 	:: ll        ! excitonic reorganization energies
        real(dp), intent(in), dimension(:,:) :: SS        ! basis transform
        complex(dpc),  dimension(:,:)        :: cc 		! correlation functions
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: temp

		real(dp), dimension(size(en,1),size(en,1)) :: S1, HH
        integer :: N, Nc, Nt
        integer :: a,b,c,m, i,j,jj, fi
        real(dp) :: kb

        integer, dimension(size(en,1),size(en,1)) :: ii
        real(dp), dimension(size(en,1),size(en,1)) :: kk,k0
        real(dp), dimension(size(en,1)) :: ga
        real(dp), dimension(size(en,1)*size(en,1)) :: oml, omr
        real(dp), dimension(size(cc,1),size(en,1)*size(en,1)) :: cw
        complex(dpc), dimension(size(cc,2),size(en,1)*size(en,1)) :: pr
        real(dp), dimension(size(cc,2)) :: x,y

        complex(dpc) :: der

		real(dp) :: rab

        real(dp) :: degen_cutoff

        character(len=128) :: char_buff

        degen_cutoff = 0.00001_dp  ! difference in energy under this value are assumed to be degenerated

        N = size(en,1)   ! number of energis
        Nc = size(cc,1)  ! number cor. fcions
        Nt = size(cc,2)  ! number of time steps in cc


        !
        ! allocate only first time you are here
        !
        if (allocated(redfield_rates_1ex)) then
            if (size(redfield_rates_1ex,1)/=N) then
                call clean_redfield()
                call allocate_redfield(N)
            end if
        else
            call allocate_redfield(N)
        end if

        if (N == 1) then

            redfield_rates_1ex(1,1) = 0.0_dp
            redfield_deph_1ex(1,1) = 0.0_dp
            redfield_deph_01ex(1) = 0.0_dp

            return
        end if


        !
        ! calculate frequencies
        !
        jj = 0;
        oml = 0
        do i = 1,N
            do j = 1,N
                 jj = jj + 1
                 ii(i,j) = jj
                 oml(jj) = en(i) - en(j)  - (ll(i) - ll(j))       ! add reorganization energies
                 omr(jj) = oml(jj)
!                 oml(jj) = abs(oml(jj))
                 if (abs(oml(jj)) <= degen_cutoff) then
                    oml(jj) = 0.0_dp
                 end if
            end do
        end do


        !
        ! calculate Fourier transforms (real parts)
        !
        do i = 1, Nt
            x(i) = (i-1)*dt
        end do
        pr = 0.0_dp

        der = 0.0_dp !(-1.0_dp,0.0_dp)/100.0_dp
        !
        ! Fourier transform of the correlation function
        !
        do i = 1, Nc

          y(1:Nt) = real(cc(i,1:Nt))
          call primitive(x,y,oml,pr)
          cw(i,:) = 2.0*real(pr(Nt,:))
          y(1:Nt) = imag(cc(i,1:Nt))
          call primitive(x,y,oml,pr)
          cw(i,:) = cw(i,:) + 2.0*real((0.0_dp,1.0_dp)*pr(Nt,:))

        end do


        !
        ! ks and gammas
        !
        do a = 1, N
            kb = 0.0_dp
            do b = 1, N
                kk(a,b) = 0.0_dp
                do m = 1, Nc
                    if (omr(ii(a,b)) >= 0.0_dp) then
                        kk(a,b) = kk(a,b) + ((abs(SS(m,a))*abs(SS(m,b)))**2)*cw(m,ii(a,b))   ! cw(m,ii(a,b) <- see May-Kuhn
                    else
                        kk(a,b) = kk(a,b) + ((abs(SS(m,a))*abs(SS(m,b)))**2)*cw(m,ii(b,a))*exp(-oml(ii(b,a))/(kB_intK*temp))
                        if (.false.) then ! if (cw(m,ii(a,b)) < 0.0_dp) then
                            write(char_buff,'(a,e12.5,a,e12.5,a,e12.5)') "Negative spectral density at  ", &
                                oml(ii(a,b)), " !  Value ", cw(m,ii(a,b)), " corrected to ", &
                                cw(m,ii(b,a))*exp(-oml(ii(b,a))/(kB_intK*temp))
                            call print_warning_message(trim(char_buff),5)
                        end if
                    end if
                 end do
                 if (a .ne. b) then
                     kb = kb + kk(a,b) ! - k0(a,b))
                 end if
            end do
            ga(a) = kb/2.0_dp !- kk(a,a)/2.0_dp
            !print *, kb
        end do

        !
        ! Construct relaxation and dephasing matrices
        !
        do a = 1, N
            do b = 1, N
                if (a == b) then
                    redfield_rates_1ex(a,b) = 0.0_dp
                    do c = 1, N
                        redfield_rates_1ex(a,a) = redfield_rates_1ex(a,a) + kk(a,c)
                    end do
                    redfield_rates_1ex(a,a) = redfield_rates_1ex(a,a) - kk(a,a)
                    redfield_deph_1ex(a,a) = 0.0_dp
                else
                    redfield_rates_1ex(a,b) = - kk(b,a)
                    redfield_deph_1ex(a,b) = ga(a) + ga(b)
                end if
            end do
            redfield_deph_01ex(a) = ga(a)
        end do


if (1 == 2) then

		!
		! saving some information
		!

        open(unit=74,file=trim(file_join("log","redfield.dat")))

        write(74,*) "# Energies: "
		do a = 1, N
			write(74,*) en(a)*Energy_internal_to_cm*Energy_cm_to_eV
		end do

        write(74,*) " "
        write(74,*) "# Redfield coeffs.:"
        write(74,*) N
        write(74,*) Nc
       	do a = 1, N
            do b = 1, N
            	rab = 0.0_dp
                do m = 1, Nc

					rab = rab + ((abs(SS(m,a))*abs(SS(m,b)))**2)

                end do
                write(74,*) rab   ! cw(m,ii(a,b) <- see May-Kuhn
            end do
         end do

        write(74, *) " "
        write(74, *) "# Basis transformation "
        do a = 1, N
          do b = 1, N
          	write(74,*) SS(a,b)
          	HH(a,b) = 0.0
          	if (a == b) then
          		HH(a,a) = en(a)
            end if
          end do
        end do

        call inv(SS,S1)
        write(74,*) " "
        do a = 1, N
          do b = 1, N
          	write(74,*) S1(a,b)
          end do
        end do

		HH = matmul(SS,matmul(HH,S1))
        write(74,*) " "
        write(74,*) "# Original Hamiltonian:"
        do a = 1, N
          do b = 1, N
          	write(74,*) HH(a,b)*Energy_internal_to_cm*Energy_cm_to_eV
          end do
        end do



        close(74)

end if


    end subroutine init_redfield_all


	!
	!  Redfield based on suplied spectral densities
	!
	subroutine init_redfield_omega(oml,SS,cw,temp)
		real(dp), intent(in), dimension(:,:) :: oml
        real(dp), intent(in), dimension(:,:) :: SS        ! basis transform
        real(dp), dimension(:,:) :: cw                     ! cw(m,ii(a,b)) = c(\omega_ab) at site m
        real(dp), intent(in) :: temp

        integer :: a,b,c,m,N,i,j,jj

        real(dp) :: kb
        integer, dimension(size(SS,1),size(SS,1)) :: ii
        real(dp), dimension(size(SS,1),size(SS,1)) :: kk,k0
        real(dp), dimension(size(SS,1)) :: ga

        real(dp) :: degen_cutoff

        character(len=128) :: char_buff

        degen_cutoff = 0.00001_dp  ! difference in energy under this value are assumed to be degenerated

        N = size(SS,1)   ! number of energis


        !
        ! allocate only first time you are here
        !
        if (allocated(redfield_rates_1ex)) then
            if (size(redfield_rates_1ex,1)/=N) then
                call clean_redfield()
                call allocate_redfield(N)
            end if
        else
            call allocate_redfield(N)
        end if

		!
		! Single transition = no rates
		!
        if (N == 1) then

            redfield_rates_1ex(1,1) = 0.0_dp
            redfield_deph_1ex(1,1) = 0.0_dp
            redfield_deph_01ex(1) = 0.0_dp
            return

        end if


        !
        ! calculate frequencies
        !
        jj = 0;
        do i = 1,N
            do j = 1,N
                 jj = jj + 1
                 ii(i,j) = jj
            end do
        end do


        !
        ! ks and gammas
        !
        do a = 1, N

            kb = 0.0_dp

            do b = 1, N
                kk(a,b) = 0.0_dp
                do m = 1, N
                    if (oml(a,b) >= 0.0_dp) then
                        kk(a,b) = kk(a,b) + ((abs(SS(m,a))*abs(SS(m,b)))**2)*cw(m,ii(a,b))   ! cw(m,ii(a,b) <- see May-Kuhn
                    else
                        kk(a,b) = kk(a,b) + ((abs(SS(m,a))*abs(SS(m,b)))**2)*cw(m,ii(b,a))*exp(-oml(b,a)/(kB_intK*temp))
                    end if
                 end do
                 if (a .ne. b) then
                     kb = kb + kk(a,b)
                 end if
            end do

            ga(a) = kb/2.0_dp

        end do

        !
        ! Construct relaxation and dephasing matrices
        !
        do a = 1, N
            do b = 1, N
                if (a == b) then
                    redfield_rates_1ex(a,b) = 0.0_dp
                    do c = 1, N
                        redfield_rates_1ex(a,a) = redfield_rates_1ex(a,a) + kk(a,c)
                    end do
                    redfield_rates_1ex(a,a) = redfield_rates_1ex(a,a) - kk(a,a)
                    redfield_deph_1ex(a,a) = 0.0_dp
                else
                    redfield_rates_1ex(a,b) = - kk(b,a)
                    redfield_deph_1ex(a,b) = ga(a) + ga(b)
                end if
            end do
            redfield_deph_01ex(a) = ga(a)
        end do

	end subroutine init_redfield_omega



    !
    ! Cleaning the module
    !
    subroutine clean_redfield()
        deallocate(redfield_rates_1ex)
        deallocate(redfield_deph_1ex)
        deallocate(redfield_deph_01ex)
    end subroutine clean_redfield

    !
    ! Allocate everything needed for calculations
    !
    subroutine allocate_redfield(N)
        integer, intent(in) :: N

        allocate(redfield_rates_1ex(N,N))
        allocate(redfield_deph_1ex(N,N))
        allocate(redfield_deph_01ex(N))

    end subroutine allocate_redfield

end module sci_redfield
