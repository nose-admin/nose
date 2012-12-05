!
!  Third order response functions
!
module response3

    use std_types
#include <util_allocation.h>
	use resp3_rss
    use util_allocation
    use resources
    use orfact
    use resp3_site
    use resp3_ext
    use resp3_wd

    implicit none

	integer		:: use_Uee_index
	logical    :: nonsecular
    character(len=512), private :: cbuff

contains


    subroutine test_response3()

        call init_response3()

        call clean_response3()

    end subroutine test_response3

    !
    !
    !
    subroutine init_response3()

        integer :: t2i
        integer(i4b) :: i, k, l, j, m, ti
		character(len=120) :: buff

!        use_Uee_index = 1 !INT(gt(2)/gt(1)) + 1

        ! writes out the time when Uee is taken
        write(buff,*) 'Uee used at T = ', (use_Uee_index-1)*dt*gt(2), ' fs, at index ', use_Uee_index
		call print_log_message(adjustl(trim(buff)),7)
		write(*,*) 'pokus'

        window_doorway = .false.

!
! **********************************************************************************
!

        ! find out how big is the problem we treat
        Ng = 0
        Ne = 0
        Nf = 0
        do i = 1, nr_blocks

            Ng = Ng + size(evops(i,i)%Ueg,2)
            Nbe = size(evops(i,i)%Ueg,1)
            Ne = Ne + Nbe

	        if(associated(evops(i,i)%Ufe)) then
    	    	Nf = size(evops(i,i)%Ufe,1)
        	elseif(associated(evops(i,i)%UfeS) ) then
	        	Nf = size(evops(i,i)%UfeS,1)
    	    else
        		Nf = 0
	        end if


        end do

        write(*,*) Ng, Ne, Nf
        write(*,*) associated(evops(1,1)%Ueg), associated(evops(1,1)%Ufe), associated(evops(1,1)%UfeS)

        ! time step
        ddt1 = dt*gt(1)

        Nt1 = size(evops(1,1)%Ueg,5) !Nt(1) !2**7
        Nt3 = Nt1

        t2i = 1

        !
        ! allocate all the quantities
        !
        ALLOCATE(rhoR,(Ng,Ne))
        ALLOCATE(rhoL,(Ne,Ng))
        ALLOCATE(dge,(Ng,Ne))
        ALLOCATE(deg,(Ne,Ng))
        ALLOCATE(dfe,(Nf,Ne))
        ALLOCATE(oeg,(Ne,Ng))
        ALLOCATE(ofe,(Nf,Ne))
        ALLOCATE(pp,(Ng))


        if (nr_blocks > 1) then
            print *, "2D is not implemented for multiple blocks"
            stop
        end if

        Np = 1
        i = 1

        Nbe = iblocks(i,i)%sblock%N1
        Nbf = Nbe*(Nbe-1)/2


        if (use_twoexcitons) then
            !print *, size(oeg,1), size(oeg,2), size(iblocks(i,i)%eblock%en)
            print *, size(deg,1), size(deg,2), size(iblocks(i,i)%eblock%dd,1), size(iblocks(i,i)%eblock%dd,2)

            do k = 1, Ne
                do l = 1, Ng

				deg(k,l) = iblocks(i,i)%eblock%dd(k,l)
                oeg(k,l) = iblocks(i,i)%eblock%en(k) - iblocks(i,i)%eblock%eng(l)

                end do
            end do

        else
        	do k = 1, Ne
                do l = 1, Ng

				deg(k,l) = iblocks(i,i)%sblock%dd(k,l)

                end do
            end do

            dfe(:,:) = 0.0_dp
            j = 0
            do l = 1, Nbe
                do k = l+1, Nbe
                    j = j + 1
                    do m = 1, Nbe
                      if (m == k) then
                        dfe(j,m) = deg(l,1)
                      else if (m == l) then
                        dfe(j,m) = deg(k,1)
                      end if
                    end do
                end do
            end do
            !dfe(Np:Np+Nbf-1,Np:Np+Nbe-1) = iblocks(i,i)%sblock%dd
            do k = 1, size(iblocks(i,i)%eblock%en,1)
                do l = 1, size(iblocks(i,i)%eblock%eng,1)

                oeg(k,l) = iblocks(i,i)%sblock%en(k) - iblocks(i,i)%eblock%eng(l)

                end do
            end do

        end if

        oeg = oeg - rwa

        ! population of the ground state
        call make_pp(pp)


        !   Np = Np + Nbe

        call make_rhoR(rhoR)
        call make_rhoL(rhoL)

        if (use_twoexcitons) then

            i = 1

            do k = 1, Nf
                do l = 1, Ne

				dfe(k,l) = iblocks(i,i)%eblock%dd_2(k,l)
                ofe(k,l) = iblocks(i,i)%eblock%en_2(k) - iblocks(i,i)%eblock%en(l)

                end do
            end do

            ofe = ofe - rwa

         end if


		! here we rename (abreviate) the evolution operators

        Uegeg => evops(1,1)%Ueg
        Ueeee => evops(1,1)%Uee
        Ugggg => evops(1,1)%Ugg
        Ufefe => evops(1,1)%Ufe


        if ( window_doorway) then

            ALLOCATE(W1g,(Ne,Ne,Nt3))
            ALLOCATE(D1g,(Ne,Ne,Nt1))
            ALLOCATE(R1g0,(Nt1,Nt3))
            call make_W1g(W1g)
            call make_D1g(D1g)
            call make_R1g(t2i)
            !DEALLOCATE(W1g)
            DEALLOCATE(D1g)

            ALLOCATE(W2g,(Ne,Ne,Nt3))
            ALLOCATE(D2g,(Ne,Ne,Nt1))
            ALLOCATE(R2g0,(Nt1,Nt3))
            W2g = W1g
            call make_D2g(D2g)
            call make_R2g(t2i)
		    DEALLOCATE(D2g)
		    DEALLOCATE(W2g)
		    DEALLOCATE(W1g)

            ALLOCATE(W3g,(Ng,Ng,Nt3))
            ALLOCATE(D3g,(Ng,Ng,Nt1))
            ALLOCATE(R3g0,(Nt1,Nt3))
            call make_W3g(W3g)
            call make_D3g(D3g)
            call make_R3g(t2i)
		    !DEALLOCATE(W3g)
		    DEALLOCATE(D3g)

            ALLOCATE(W4g,(Ne,Ne,Nt3))
            ALLOCATE(D4g,(Ne,Ne,Nt1))
            ALLOCATE(R4g0,(Nt1,Nt3))
            W4g = W3g
            call make_D4g(D4g)
            call make_R4g(t2i)
		    DEALLOCATE(W4g)
		    DEALLOCATE(D4g)
		    DEALLOCATE(W3g)

            ALLOCATE(W1f,(Ne,Ne,Nt3))
            ALLOCATE(D1f,(Ne,Ne,Nt1))
            ALLOCATE(R1f0,(Nt1,Nt3))
            call make_W1f(W1f)
            call make_D1f(D1f)
            call make_R1f(t2i)
		    !DEALLOCATE(W1f)
		    DEALLOCATE(D1f)

	        ALLOCATE(W2f,(Ne,Ne,Nt3))
    	    ALLOCATE(D2f,(Ne,Ne,Nt1))
        	ALLOCATE(R2f0,(Nt1,Nt3))
        	W2f = W1f
        	call make_D2f(D2f)
        	call make_R2f(t2i)
			DEALLOCATE(D2f)
			DEALLOCATE(W2f)
			DEALLOCATE(W1f)


        else

			!nonsecular = .false.

			if (nonsecular) then
			print *, "R1g "

            !
            ! R1g
            !
            ALLOCATE(extD,(Ne,Ng,Ne,Ng,Nt1))
            ALLOCATE(R1g0,(Nt1,Nt3))
            call make_orD1g(extD)
            call make_orR1g(t2i)
!            DEALLOCATE(extD)

			print *, "R2g "
			!
            ! R2g
            !
!            ALLOCATE(extD,(Ne,Ng,Ne,Ng,Nt1))
            ALLOCATE(R2g0,(Nt1,Nt3))
            call make_orD2g(extD)
            call make_orR2g(t2i)
!            DEALLOCATE(extD)

			print *, "R3g "
            !
            ! R3g
            !
!            ALLOCATE(extD,(Ne,Ng,Ne,Ng,Nt1))
            ALLOCATE(R3g0,(Nt1,Nt3))
            call make_orD3g(extD)
            call make_orR3g(t2i)
!            DEALLOCATE(extD)

			print *, "R4g "
            !
            ! R4g
            !
!            ALLOCATE(extD,(Ne,Ng,Ne,Ng,Nt1))
            ALLOCATE(R4g0,(Nt1,Nt3))
            call make_orD4g(extD)
            call make_orR4g(t2i)
            DEALLOCATE(extD)

			print *, "R1f "

			!
			! R1f
			!
            ALLOCATE(extD,(Ne,Nf,Ne,Nf,Nt1))
            ALLOCATE(R1f0,(Nt1,Nt3))
            print *, "... D1f"
            call make_orD1f(extD)
            print *, "... R1f"
            call make_orR1f(t2i)
!            DEALLOCATE(extD)

			print *, "R2f "
			!
			! R2f
			!
!            ALLOCATE(extD,(Ne,Nf,Ne,Nf,Nt1))
            ALLOCATE(R2f0,(Nt1,Nt3))
            print *, "... D2f"
            call make_orD2f(extD)
            print *, "... R2f"
            call make_orR2f(t2i)
            DEALLOCATE(extD)

			else !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			write(cbuff,'(a,2f18.12)') "R1g secular "
            call print_log_message(trim(cbuff),7)

            !
            ! R1g
            !
            print *, "R1g secular "
            ALLOCATE(extD,(Ne,Ng,Ne,Ng,Nt1))
            ALLOCATE(R1g0,(Nt1,Nt3))
            call make_orD1g_sec(extD)
            call make_orR1g_sec(t2i)
!            DEALLOCATE(extD)

			write(cbuff,'(a,2f18.12)') "R2g secular "
            call print_log_message(trim(cbuff),7)
			!
            ! R2g
            !
            print *, "R2g secular "
!            ALLOCATE(extD,(Ne,Ng,Ne,Ng,Nt1))
            ALLOCATE(R2g0,(Nt1,Nt3))
            call make_orD2g_sec(extD)
            call make_orR2g_sec(t2i)
!            DEALLOCATE(extD)

			write(cbuff,'(a,2f18.12)') "R3g secular "
            call print_log_message(trim(cbuff),7)
            !
            ! R3g
            !
            print *, "R3g secular "
!            ALLOCATE(extD,(Ne,Ng,Ne,Ng,Nt1))
            ALLOCATE(R3g0,(Nt1,Nt3))
            call make_orD3g_sec(extD)
            call make_orR3g_sec(t2i)
!            DEALLOCATE(extD)

			write(cbuff,'(a,2f18.12)') "R4g secular "
            call print_log_message(trim(cbuff),7)
            !
            ! R4g
            !
            print *, "R4g secular "
!            ALLOCATE(extD,(Ne,Ng,Ne,Ng,Nt1))
            ALLOCATE(R4g0,(Nt1,Nt3))
            call make_orD4g_sec(extD)
            call make_orR4g_sec(t2i)
            DEALLOCATE(extD)

			write(cbuff,'(a,2f18.12)') "R1f secular "
            call print_log_message(trim(cbuff),7)

			!
			! R1f
			!
            ALLOCATE(extD,(Ne,Nf,Ne,Nf,1))
            ALLOCATE(R1f0,(Nt1,Nt3))
!            ALLOCATE(extD,(Ne,Nf,Ne,Nf,Nt1))
            ALLOCATE(R2f0,(Nt1,Nt3))

            do ti = 1, Nt1 ! significantly saves memory
            	call make_orD1f_sec(extD,ti)
            	call make_orR1f_sec(t2i, ti)

				write(cbuff,'(a,2f18.12)') "R2f secular "
            	call print_log_message(trim(cbuff),7)
            !
			! R2f
			!
            	call make_orD2f_sec(extD,ti)
            	call make_orR2f_sec(t2i,ti)
            end do

            DEALLOCATE(extD)


			end if

        end if



if (.false.) then

        if ( window_doorway) then


        ALLOCATE(W3f,(Ng,Nf,Nt3))
        ALLOCATE(D3f,(Ng,Nf,Nt1))
        ALLOCATE(R3f0,(Nt1,Nt3))
        call make_W3f(W3f)
        call make_D3f(D3f)
        call make_R3f(t2i)
		DEALLOCATE(D3f)
		DEALLOCATE(W3f)

        else

        ALLOCATE(W3f,(Ng,Nf,Nt3))
        ALLOCATE(D3f,(Ng,Nf,Nt1))
        ALLOCATE(R3f0,(Nt1,Nt3))
        call make_W3f(W3f)
        call make_D3f(D3f)
        call make_R3f(t2i)
        DEALLOCATE(D3f)
        DEALLOCATE(W3f)

        end if


        if ( window_doorway) then


        ALLOCATE(W4f,(Nf,Ng,Nt3))
        ALLOCATE(D4f,(Nf,Ng,Nt1))
        ALLOCATE(R4f0,(Nt1,Nt3))
        call make_W4f(W4f)
        call make_D4f(D4f)
        call make_R4f(t2i)
		DEALLOCATE(D4f)
		DEALLOCATE(W4f)

        else

        ALLOCATE(W4f,(Nf,Ng,Nt3))
        ALLOCATE(D4f,(Nf,Ng,Nt1))
        ALLOCATE(R4f0,(Nt1,Nt3))
        call make_W4f(W4f)
        call make_D4f(D4f)
        call make_R4f(t2i)
        DEALLOCATE(D4f)
        DEALLOCATE(W4f)

        end if


end if

    end subroutine init_response3

    !
    !
    !
    subroutine clean_response3()

        DEALLOCATE(rhoR)
        DEALLOCATE(rhoL)
        DEALLOCATE(dge)
        DEALLOCATE(deg)
        DEALLOCATE(dfe)
        DEALLOCATE(pp)
        DEALLOCATE(oeg)
        DEALLOCATE(ofe)


        DEALLOCATE(R1g0)
        DEALLOCATE(R2g0)
        DEALLOCATE(R3g0)
        DEALLOCATE(R4g0)
        DEALLOCATE(R1f0)
        DEALLOCATE(R2f0)
        !DEALLOCATE(R3f0)
        !DEALLOCATE(R4f0)



    end subroutine clean_response3


    subroutine make_rhoL(r)
        real(dp), dimension(:,:), intent(out) :: r
        ! local
        integer :: i,j

        do i = 1,Ne
            do j = 1, Ng
                rhoL(i,j) = deg(i,j)*pp(j)
            end do
        end do

    end subroutine make_rhoL

    subroutine make_rhoR(r)
        real(dp), dimension(:,:), intent(out) :: r
        ! local
        integer :: i,j

        do i = 1,Ng
            do j = 1, Ne
                rhoR(i,j) = pp(i)*deg(j,i)
            end do
        end do

    end subroutine make_rhoR

    !
    ! Equilibrium population of the ground state
    !
    subroutine make_pp(p)
        real(dp), dimension(:), intent(out) :: p
        ! local
        integer :: i, N
        real(dp) :: sum

        if(size(p) /= size(evops(1,1)%Ugg,1)) then
        	call print_error_message(-1, 'dimension error in make_pp, module response3')
        end if

        N = size(evops(1,1)%Ugg,1)

        if (N == 1) then
            p(1) = 1.0_dp
        else
            sum = 0.0_dp
            do i = 1, N
                p(i) = exp(-iblocks(1,1)%eblock%eng(i) / temp / kB_intK)
                sum = sum + p(i)
            end do
            p = p/sum
        end if

    end subroutine make_pp




    !
    ! orD1g propagator
    !
    subroutine make_orD1g(DD)
        complex(dpc), dimension(:,:,:,:,:), intent(out) :: DD
        ! local
        integer :: c,e,d,b,a,f,g,h,i,j,ti

        DD = 0.0_dp
        do ti = 1, Nt1

            do i = 1, Ne
            do j = 1, Ng
            do f = 1, Ne
            do h = 1, Ng

                do g = 1, Ne
                do d = 1, Ng
                do e = 1, Ne
                do c = 1, Ne
                do b = 1, Ne
                do a = 1, Ng
                    DD(i,j,f,h,ti) = DD(i,j,f,h,ti) + &          ! def(g,h)
                     orfact_23_S(i,j,g,h,e,d,b,a)* &
                    deg(i,j)*deg(g,h)*Ueeee(f,g,c,e,use_Uee_index) &
                    *Uegeg(c,d,b,a,ti)*rhoL(b,a)*deg(e,d)
                                                ! dge(d,e)

                end do
                end do
                end do
                end do
                end do
                end do

            end do
            end do
            end do
            end do

        end do

    end subroutine make_orD1g


    !
    ! orD1g propagator  SECULAR VERSION
    !
    subroutine make_orD1g_sec(DD)
        complex(dpc), dimension(:,:,:,:,:), intent(out) :: DD
        ! local
        integer :: c,e,d,b,a,f,g,h,i,j,ti

        DD = 0.0_dp
        do ti = 1, Nt1

			! populations
            do i = 1, Ne
            do j = 1, Ng
            !do f = 1, Ne
            f = i
            do h = 1, Ng

                !do g = 1, Ne
                g = f
                do d = 1, Ng
                !do e = 1, Ne
                !do c = 1, Ne
                do b = 1, Ne
                c = b
                e = c
                do a = 1, Ng
                    DD(i,j,f,h,ti) = DD(i,j,f,h,ti) + &          ! def(g,h)
                     orfact_23_S(i,j,g,h,e,d,b,a)* &
                    deg(i,j)*deg(g,h)*Ueeee(f,g,c,e,use_Uee_index) &
                    *Uegeg(c,d,b,a,ti)*rhoL(b,a)*deg(e,d)
                                                ! dge(d,e)

                end do
                end do
                !end do
                !end do
                end do
                !end do

            end do
            !end do
            end do
            end do

			! coherences
            do i = 1, Ne
            do j = 1, Ng
            !do f = 1, Ne
            f = i
            do h = 1, Ng

                do g = 1, Ne
                if (f /= g) then
                do d = 1, Ng
                !do e = 1, Ne
                !do c = 1, Ne
                e = g
                c = f
                !do b = 1, Ne
                b = c
                do a = 1, Ng
                    DD(i,j,f,h,ti) = DD(i,j,f,h,ti) + &          ! def(g,h)
                     orfact_23_S(i,j,g,h,e,d,b,a)* &
                    deg(i,j)*deg(g,h)*Ueeee(f,g,c,e,use_Uee_index) &
                    *Uegeg(c,d,b,a,ti)*rhoL(b,a)*deg(e,d)
                                                ! dge(d,e)

                end do
                !end do
                !end do
                !end do
                end do
                end if
                end do

            end do
            !end do
            end do
            end do


        end do

    end subroutine make_orD1g_sec


    !
    ! orD1f propagator
    !
    subroutine make_orD1f(DD)
        complex(dpc), dimension(:,:,:,:,:), intent(out) :: DD
        ! local
        integer :: c,e,d,b,a,f,g,h,i,j,ti

        DD = 0.0_dp
        do ti = 1, Nt1

            do i = 1, Ne
            do j = 1, Nf
            do f = 1, Ne
            !do h = 1, Nf
            h = j

                do g = 1, Ne
                do d = 1, Ng
                do e = 1, Ne
                do c = 1, Ne
             !   do b = 1, Ne
                b = c
                do a = 1, Ng

                    DD(i,j,f,h,ti) = DD(i,j,f,h,ti) + &          ! def(g,h)
                     orfact_e_S(j,i,h,g,e,d,b,a)* &
                    dfe(j,i)*dfe(h,g)*Ueeee(f,g,c,e,use_Uee_index) &
                    *Uegeg(c,d,b,a,ti)*rhoL(b,a)*deg(e,d)
                                                ! dge(d,e)
                end do
                !end do
                end do
                end do
                end do
                end do

            !end do
            end do
            end do
            end do

        end do

    end subroutine make_orD1f


    !
    ! orD1f propagator   SECULAR VERSION
    !
    subroutine make_orD1f_sec(DD, ti)
        complex(dpc), dimension(:,:,:,:,:), intent(out) :: DD
        integer(i4b), intent(in) :: ti
        ! local
        integer(i4b) :: c,e,d,b,a,f,g,h,i,j

        DD = 0.0_dp
        !do ti = 1, Nt1

		    ! coherences

            do i = 1, Ne
            do j = 1, Nf
            !do f = 1, Ne
            f = i
            !do h = 1, Nf
            h = j


                do g = 1, Ne
                if (f/=g) then
                do d = 1, Ng
             !   do e = 1, Ne
             !   do c = 1, Ne
             !   do b = 1, Ne
                e = g
                c = f
                b = c
                a = d
                !do a = 1, Ng

!                    DD(i,j,f,h,ti) = DD(i,j,f,h,ti) + &          ! def(g,h)
                    DD(i,j,f,h,1) = DD(i,j,f,h,1) + &          ! def(g,h)
                     orfact_e_S(j,i,h,g,e,d,b,a)* &
                    dfe(j,i)*dfe(h,g)*Ueeee(f,g,c,e,use_Uee_index) &
                    *Uegeg(c,d,b,a,ti)*rhoL(b,a)*deg(e,d)
                                                ! dge(d,e)
                !end do
                !end do
                !end do
                !end do
                end do
                end if
                end do


            !end do
            !end do
            end do
            end do

		    ! populations

            do i = 1, Ne
            do j = 1, Nf
            !do f = 1, Ne
            f = i
            !do h = 1, Nf
            h = j


             !   do g = 1, Ne
                g = f
                do d = 1, Ng
                do e = 1, Ne
             !   do c = 1, Ne
             !   do b = 1, Ne
                c = e
                b = c
                a = d
                !do a = 1, Ng

!                    DD(i,j,f,h,ti) = DD(i,j,f,h,ti) + &          ! def(g,h)
                    DD(i,j,f,h,1) = DD(i,j,f,h,1) + &          ! def(g,h)
                     orfact_e_S(j,i,h,g,e,d,b,a)* &
                    dfe(j,i)*dfe(h,g)*Ueeee(f,g,c,e,use_Uee_index) &
                    *Uegeg(c,d,b,a,ti)*rhoL(b,a)*deg(e,d)
                                                ! dge(d,e)
                !end do
                !end do
                !end do
                end do
                end do
                !end do


            !end do
            !end do
            end do
            end do


        !end do

    end subroutine make_orD1f_sec


    !
    ! orD2g propagator
    !
    subroutine make_orD2g(DD)
        complex(dpc), dimension(:,:,:,:,:), intent(out) :: DD
        ! local
        integer :: c,e,d,b,a,f,g,h,i,j,ti

        DD = 0.0_dp
        do ti = 1, Nt1

            do i = 1, Ne
            do j = 1, Ng
            do f = 1, Ne
            do h = 1, Ng

                do g = 1, Ne
                do d = 1, Ne
                do e = 1, Ne
                do c = 1, Ng
                do b = 1, Ne
                do a = 1, Ng
                    DD(i,j,f,h,ti) = DD(i,j,f,h,ti) + &
                        orfact_23_S(i,j,g,h,e,c,b,a)* &
                       	deg(i,j)*Ueeee(f,g,e,d,use_Uee_index)*deg(e,c) &
                       	*conjg(Uegeg(d,c,b,a,ti))*rhoR(a,b)*deg(g,h) !Ugege(c,d,a,b,ti)*rhoR(a,b)

                end do
                end do
                end do
                end do
                end do
                end do

                !stop

            end do
            end do
            end do
            end do

        end do



    end subroutine make_orD2g


    !
    ! orD2g propagator  SECULAR VERSION
    !
    subroutine make_orD2g_sec(DD)
        complex(dpc), dimension(:,:,:,:,:), intent(out) :: DD
        ! local
        integer :: c,e,d,b,a,f,g,h,i,j,ti

		print *,"make_orD2g_sec"

        DD = 0.0_dp
        do ti = 1, Nt1

			! coherences
            do i = 1, Ne
            do j = 1, Ng
            do f = 1, Ne
            do h = 1, Ng

                do g = 1, Ne
                if (f /= g) then
                !do d = 1, Ne
                !do e = 1, Ne
                e = f
                d = g
                do c = 1, Ng
                !do b = 1, Ne
                b = d
                !do a = 1, Ng
                a = c
                    DD(i,j,f,h,ti) = DD(i,j,f,h,ti) + &
                        orfact_23_S(i,j,g,h,e,c,b,a)* &
                       	deg(i,j)*Ueeee(f,g,e,d,use_Uee_index)*deg(e,c) &
                       	*conjg(Uegeg(d,c,b,a,ti))*rhoR(a,b)*deg(g,h) !Ugege(c,d,a,b,ti)*rhoR(a,b)

                !end do
                !end do
                !end do
                !end do
                end do
                end if
                end do

            end do
            end do
            end do
            end do

			! populations
            do i = 1, Ne
            do j = 1, Ng
            do f = 1, Ne
            do h = 1, Ng

                !do g = 1, Ne
                do d = 1, Ne
                !do e = 1, Ne
                g = f
                e = d
                do c = 1, Ng
                !do b = 1, Ne
                b = d
                !do a = 1, Ng
                a = c
                    DD(i,j,f,h,ti) = DD(i,j,f,h,ti) + &
                        orfact_23_S(i,j,g,h,e,c,b,a)* &
                       	deg(i,j)*Ueeee(f,g,e,d,use_Uee_index)*deg(e,c) &
                       	*conjg(Uegeg(d,c,b,a,ti))*rhoR(a,b)*deg(g,h) !Ugege(c,d,a,b,ti)*rhoR(a,b)

                !end do
                !end do
                !end do
                !end do
                end do
                end do

            end do
            end do
            end do
            end do


        end do



    end subroutine make_orD2g_sec


    !
    ! D2f propagator
    !
    !subroutine make_orD2f(DD)
    !    complex(dpc), dimension(:,:,:), intent(out) :: DD
    !
    !    call make_D2g(DD)
    !
    !end subroutine make_orD2f
    subroutine make_orD2f(DD)
        complex(dpc), dimension(:,:,:,:,:), intent(out) :: DD
        ! local
        integer :: c,e,d,b,a,f,g,h,i,j,ti

        DD = 0.0_dp
        do ti = 1, Nt1

            do i = 1, Ne
            do j = 1, Nf
            do f = 1, Ne
            !do h = 1, Nf
            h = j

                do g = 1, Ne
                do e = 1, Ne
                do d = 1, Ne
                do c = 1, Ng
                b = d
                !do b = 1, Ne
                do a = 1, Ng
                    DD(i,j,f,h,ti) = DD(i,j,f,h,ti) + &          ! def(g,h)
                     orfact_e_S(j,i,h,g,e,c,b,a)* &
                    dfe(j,i)*dfe(h,g)*Ueeee(f,g,e,d,use_Uee_index) &
                    *conjg(Uegeg(d,c,b,a,ti))*rhoR(a,b)*deg(e,c)
                                                ! dge(d,e)
                end do
                !end do
                end do
                end do
                end do
                end do

            !end do
            end do
            end do
            end do

        end do

    end subroutine make_orD2f


    !
    ! D2f propagator  SECULAR VERSION
    !
    !subroutine make_orD2f(DD)
    !    complex(dpc), dimension(:,:,:), intent(out) :: DD
    !
    !    call make_D2g(DD)
    !
    !end subroutine make_orD2f
    subroutine make_orD2f_sec(DD,ti)
        complex(dpc), dimension(:,:,:,:,:), intent(out) :: DD
        integer(i4b), intent(in) :: ti
        ! local
        integer(i4b) :: c,e,d,b,a,f,g,h,i,j

        DD = 0.0_dp
        !do ti = 1, Nt1

			! populations
            do i = 1, Ne
            do j = 1, Nf
            !do f = 1, Ne
            f = i
            !do h = 1, Nf
            h = j

                !do g = 1, Ne
                do e = 1, Ne
                !do d = 1, Ne
                d = e
                do c = 1, Ng
                g = f
                b = d
                !do b = 1, Ne
                !do a = 1, Ng
                a = c
!                    DD(i,j,f,h,ti) = DD(i,j,f,h,ti) + &          ! def(g,h)
                    DD(i,j,f,h,1) = DD(i,j,f,h,1) + &          ! def(g,h)
                     orfact_e_S(j,i,h,g,e,c,b,a)* &
                    dfe(j,i)*dfe(h,g)*Ueeee(f,g,e,d,use_Uee_index) &
                    *conjg(Uegeg(d,c,b,a,ti))*rhoR(a,b)*deg(e,c)
                                                ! dge(d,e)
                !end do
                !end do
                end do
                !end do
                end do
                !end do

            !end do
            !end do
            end do
            end do

			! coherence
            do i = 1, Ne
            do j = 1, Nf
            !do f = 1, Ne
            f = i
            !do h = 1, Nf
            h = j

                do g = 1, Ne
                if (f /= g) then
                !do e = 1, Ne
                e = f
                !do d = 1, Ne
                d = g
                do c = 1, Ng
                b = d
                !do b = 1, Ne
                !do a = 1, Ng
                a = c
!                    DD(i,j,f,h,ti) = DD(i,j,f,h,ti) + &          ! def(g,h)
                    DD(i,j,f,h,1) = DD(i,j,f,h,1) + &          ! def(g,h)
                     orfact_e_S(j,i,h,g,e,c,b,a)* &
                    dfe(j,i)*dfe(h,g)*Ueeee(f,g,e,d,use_Uee_index) &
                    *conjg(Uegeg(d,c,b,a,ti))*rhoR(a,b)*deg(e,c)
                                                ! dge(d,e)
                !end do
                !end do
                end do
                !end do
                !end do
                end if
                end do

            !end do
            !end do
            end do
            end do

        !end do

    end subroutine make_orD2f_sec


    !
    ! D3g propagator
    !
    subroutine make_orD3g(DD)
        complex(dpc), dimension(:,:,:,:,:), intent(out) :: DD
        ! local
        integer :: c,e,d,b,a,f,g,h,i,j,ti

        DD = 0.0_dp
        do ti = 1, Nt1

            do i = 1, Ne
            do j = 1, Ng
            do h = 1, Ne
            do g = 1, Ng

                do f = 1, Ng
                do d = 1, Ne
                do e = 1, Ng
                do c = 1, Ng
                do b = 1, Ne
                do a = 1, Ng
                    !DD(i,j,f,h,ti) = DD(i,j,f,h,ti) &
                    !   + orfact_23(i,1,g,1,e,1,b,1)*deg(i,j)*Ueeee(f,g,e,d,use_Uee_index)*deg(e,c) &
                    !   *conjg(Uegeg(d,c,b,a,ti))*rhoR(a,b)*deg(g,h) !Ugege(c,d,a,b,ti)*rhoR(a,b)
                    !print *, orfact_23(i,1,g,1,e,1,b,1)

                    DD(i,j,h,g,ti) = DD(i,j,h,g,ti) + &
                        orfact_23_S(i,j,h,f,d,e,b,a)* &
                        deg(i,j)*deg(h,f)*Ugggg(f,g,c,e,use_Uee_index)*  &
                        conjg(Uegeg(d,c,b,a,ti))*rhoR(a,b)*deg(d,e)
                end do
                end do
                end do
                end do
                end do
                end do

                !stop

            end do
            end do
            end do
            end do

        end do


    end subroutine make_orD3g


    !
    ! D3g propagator  SECULAR VERSION
    !
    subroutine make_orD3g_sec(DD)
        complex(dpc), dimension(:,:,:,:,:), intent(out) :: DD
        ! local
        integer :: c,e,d,b,a,f,g,h,i,j,ti

        DD = 0.0_dp
        do ti = 1, Nt1

            do i = 1, Ne
            do j = 1, Ng
            !do h = 1, Ne
            h = i
            !do g = 1, Ng
			g = j

                do f = 1, Ng
                do d = 1, Ne
                do e = 1, Ng
                do c = 1, Ng
                !do b = 1, Ne
                b = d
                !do a = 1, Ng
                a = c
                    DD(i,j,h,g,ti) = DD(i,j,h,g,ti) + &
                        orfact_23_S(i,j,h,f,d,e,b,a)* &
                        deg(i,j)*deg(h,f)*Ugggg(f,g,c,e,use_Uee_index)*  &
                        conjg(Uegeg(d,c,b,a,ti))*rhoR(a,b)*deg(d,e)
                !end do
                !end do
                end do
                end do
                end do
                end do

                !stop

            !end do
            !end do
            end do
            end do

        end do


    end subroutine make_orD3g_sec



    !
    ! D3g propagator
    !
    subroutine make_orD4g(DD)
        complex(dpc), dimension(:,:,:,:,:), intent(out) :: DD
        ! local
        integer :: c,e,d,b,a,f,g,h,i,j,ti

        DD = 0.0_dp
        do ti = 1, Nt1

            do i = 1, Ne
            do j = 1, Ng
            do h = 1, Ne
            do g = 1, Ng

                do f = 1, Ng
                do d = 1, Ng
                do e = 1, Ng
                do c = 1, Ne
                do b = 1, Ne
                do a = 1, Ng

                    DD(i,j,h,g,ti) = DD(i,j,h,g,ti) + &
                        orfact_23_S(i,j,h,f,c,e,b,a)* &
                        deg(i,j)*deg(h,f)*Ugggg(f,g,e,d,use_Uee_index)*  &
                        Uegeg(c,d,b,a,ti)*rhoL(b,a)*deg(c,e)
                end do
                end do
                end do
                end do
                end do
                end do

                !stop

            end do
            end do
            end do
            end do

        end do


    end subroutine make_orD4g

    !
    ! D3g propagator   SECULAR VERSION
    !
    subroutine make_orD4g_sec(DD)
        complex(dpc), dimension(:,:,:,:,:), intent(out) :: DD
        ! local
        integer :: c,e,d,b,a,f,g,h,i,j,ti

        DD = 0.0_dp
        do ti = 1, Nt1

			! coherences
            do i = 1, Ne
            do j = 1, Ng
            do h = 1, Ne
            do g = 1, Ng

                do f = 1, Ng
                if (f /= g) then
                !do d = 1, Ng
                d = g
                !do e = 1, Ng
                e = f
                do c = 1, Ne
                !do b = 1, Ne
                b = c
                !do a = 1, Ng
                a = d

                    DD(i,j,h,g,ti) = DD(i,j,h,g,ti) + &
                        orfact_23_S(i,j,h,f,c,e,b,a)* &
                        deg(i,j)*deg(h,f)*Ugggg(f,g,e,d,use_Uee_index)*  &
                        Uegeg(c,d,b,a,ti)*rhoL(b,a)*deg(c,e)
                !end do
                !end do
                !end do
                !end do
                end do
                end if
                end do

                !stop

            end do
            end do
            end do
            end do

			! populations
            do i = 1, Ne
            do j = 1, Ng
            do h = 1, Ne
            do g = 1, Ng

                !do f = 1, Ng
                f = g
                do d = 1, Ng
                !do e = 1, Ng
                e = d
                do c = 1, Ne
                !do b = 1, Ne
                b = c
                !do a = 1, Ng
                a = d

                    DD(i,j,h,g,ti) = DD(i,j,h,g,ti) + &
                        orfact_23_S(i,j,h,f,c,e,b,a)* &
                        deg(i,j)*deg(h,f)*Ugggg(f,g,e,d,use_Uee_index)*  &
                        Uegeg(c,d,b,a,ti)*rhoL(b,a)*deg(c,e)
                !end do
                !end do
                !end do
                !end do
                end do
                end do

                !stop

            end do
            end do
            end do
            end do

        end do


    end subroutine make_orD4g_sec



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    ! Liouville Pathways
    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  R1g
    !
    subroutine make_R1g(ti)
        integer :: ti
        ! local
        integer :: t1, t3, a, b

        if (ti == 1) then

            do t1 = 1, Nt1
                do t3 = 1, Nt1
                    R1g0(t3,t1) = 0.0_dp
                    do a = 1, Ne
                        do b = 1, Ne
                            R1g0(t3,t1) = R1g0(t3,t1) + W1g(a,b,t3)*D1g(a,b,t1)
                        end do
                    end do
                end do
            end do

        else

        end if


    end subroutine make_R1g

    !
    !     R1g with orientational averaging
    !
    subroutine make_orR1g(ti)
        integer :: ti

        ! local
        integer :: t1, t3, i,j,f,h

        if (ti == 1) then

            R1g0 = 0.0_dp

            do t1 = 1, Nt1
                do t3 = 1, Nt1

                    do i = 1, Ne
                        do j = 1, Ng
                            do f = 1, Ne
                                do h = 1, Ng
                                    R1g0(t3,t1) = R1g0(t3,t1) + evops(1,1)%Ueg(i,j,f,h,t3)*extD(i,j,f,h,t1)
                                                               !evops(1,1)%Ufe(i,j,f,h,t3)*extD(i,j,f,h,t1) !W2g(a,b,t3)*D2g(a,b,t1)
                                end do
                            end do
                        end do
                    end do

                end do
            end do

        else

        end if

    end subroutine make_orR1g

    !
    !     R1g with orientational averaging  SECULAR VERSION
    !
    subroutine make_orR1g_sec(ti)
        integer :: ti

        ! local
        integer :: t1, t3, i,j,f,h

        if (ti == 1) then

            R1g0 = 0.0_dp

            do t1 = 1, Nt1
                do t3 = 1, Nt1

                    do i = 1, Ne
                        do j = 1, Ng
                            !do f = 1, Ne
                            f = i
                                do h = 1, Ng
                                    R1g0(t3,t1) = R1g0(t3,t1) + evops(1,1)%Ueg(i,j,f,h,t3)*extD(i,j,f,h,t1)
                                                               !evops(1,1)%Ufe(i,j,f,h,t3)*extD(i,j,f,h,t1) !W2g(a,b,t3)*D2g(a,b,t1)
                                end do
                            !end do
                        end do
                    end do

                end do
            end do

        else

        end if

    end subroutine make_orR1g_sec


    !
    ! R1f
    !
    subroutine make_R1f(ti)
        integer :: ti
        ! local
        integer :: t1, t3, a, b

        if (ti == 1) then

            do t1 = 1, Nt1
                do t3 = 1, Nt1
                    R1f0(t3,t1) = 0.0_dp
                    do a = 1, Ne
                        do b = 1, Ne
                            R1f0(t3,t1) = R1f0(t3,t1) + W1f(a,b,t3)*D1f(a,b,t1)
                        end do
                    end do
                end do
            end do

        else

        end if

    end subroutine make_R1f

    !
    ! R1f with orientational averaging
    !
    subroutine make_orR1f(ti)
        integer :: ti

        ! local
        integer :: t1, t3, i,j,f,h

        if (ti == 1) then

            R1f0 = 0.0_dp

            do t1 = 1, Nt1
                do t3 = 1, Nt1

                    do i = 1, Ne
                        do j = 1, Nf
                            do f = 1, Ne
                                h = j
                                !do h = 1, Nf
                                    R1f0(t3,t1) = R1f0(t3,t1) + evops(1,1)%Ufe(j,i,h,f,t3)*extD(i,j,f,h,t1)
                                                               !evops(1,1)%Ufe(i,j,f,h,t3)*extD(i,j,f,h,t1) !W2g(a,b,t3)*D2g(a,b,t1)
                                !end do
                            end do
                        end do
                    end do

                end do
            end do

        else

        end if

    end subroutine make_orR1f


    !
    ! R1f with orientational averaging  SECULAR VERSION
    !
    subroutine make_orR1f_sec(ti,t1)
        integer(i4b), intent(in) :: ti,t1

        ! local
        integer :: t3, i,j,f,h

        if (ti == 1) then

            !R1f0 = 0.0_dp
            R1f0(:,t1) = 0.0_dp
            !write(*,*) 'make_orR1f_sec'

            !do t1 = 1, Nt1
                do t3 = 1, Nt1

                    do i = 1, Ne
                        do j = 1, Nf
                            !do f = 1, Ne
                                h = j
                                f = i
                                !do h = 1, Nf
!                                    R1f0(t3,t1) = R1f0(t3,t1) + evops(1,1)%UfeS(j,i,t3)*extD(i,j,f,h,t1)
                                    R1f0(t3,t1) = R1f0(t3,t1) + evops(1,1)%UfeS(j,i,t3)*extD(i,j,f,h,1)
                                    							! evops(1,1)%Ufe(j,i,h,f,t3)*extD(i,j,f,h,t1)
                                                               !evops(1,1)%Ufe(i,j,f,h,t3)*extD(i,j,f,h,t1) !W2g(a,b,t3)*D2g(a,b,t1)
                                !end do
                            !end do
                        end do
                    end do

                end do
            !end do

        else

        end if

    end subroutine make_orR1f_sec


    !
    ! R2g
    !
    subroutine make_R2g(ti)
        integer :: ti

        ! local
        integer :: t1, t3, a, b

        if (ti == 1) then

            do t1 = 1, Nt1
                do t3 = 1, Nt1
                    R2g0(t3,t1) = 0.0_dp
                    do a = 1, Ne
                        do b = 1, Ne
                            R2g0(t3,t1) = R2g0(t3,t1) + W2g(a,b,t3)*D2g(a,b,t1)
                        end do
                    end do
                end do
            end do

        else

        end if

    end subroutine make_R2g


    !
    ! R2g with orientational averaging
    !
    subroutine make_orR2g(ti)
        integer :: ti

        ! local
        integer :: t1, t3, i,j,f,h

        if (ti == 1) then

            R2g0 = 0.0_dp

            do t1 = 1, Nt1
                do t3 = 1, Nt1

                    do i = 1, Ne
                        do j = 1, Ng
                            do f = 1, Ne
                                do h = 1, Ng
                                    R2g0(t3,t1) = R2g0(t3,t1) + evops(1,1)%Ueg(i,j,f,h,t3)*extD(i,j,f,h,t1) !W2g(a,b,t3)*D2g(a,b,t1)
                                end do
                            end do
                        end do
                    end do

                end do
            end do

        else

        end if

    end subroutine make_orR2g


    !
    ! R2g with orientational averaging SECULAR VERSION
    !
    subroutine make_orR2g_sec(ti)
        integer :: ti

        ! local
        integer :: t1, t3, i,j,f,h

        if (ti == 1) then

            R2g0 = 0.0_dp

            do t1 = 1, Nt1
                do t3 = 1, Nt1

                    do i = 1, Ne
                        do j = 1, Ng
                            !do f = 1, Ne
                            f = i
                                do h = 1, Ng
                                    R2g0(t3,t1) = R2g0(t3,t1) + evops(1,1)%Ueg(i,j,f,h,t3)*extD(i,j,f,h,t1) !W2g(a,b,t3)*D2g(a,b,t1)
                                end do
                            !end do
                        end do
                    end do

                end do
            end do

        else

        end if

    end subroutine make_orR2g_sec


    !
    ! R2f
    !
    subroutine make_R2f(ti)
        integer :: ti
        ! local
        integer :: t1, t3, a, b

        if (ti == 1) then

            do t1 = 1, Nt1
                do t3 = 1, Nt1
                    R2f0(t3,t1) = 0.0_dp
                    do a = 1, Ne
                        do b = 1, Ne
                            R2f0(t3,t1) = R2f0(t3,t1) + W2f(a,b,t3)*D2f(a,b,t1)
                        end do
                    end do
                end do
            end do

        else

        end if

    end subroutine make_R2f

    !
    ! R2f with orientational averaging
    !
    subroutine make_orR2f(ti)
        integer :: ti

        ! local
        integer :: t1, t3, i,j,f,h

        if (ti == 1) then

            R2f0 = 0.0_dp

            do t1 = 1, Nt1
                do t3 = 1, Nt1

                    do i = 1,Ne
                        do j = 1, Nf
                            do f = 1, Ne
                                !do h = 1, Nf
                                 !h = 3
                                 h = j
                                    R2f0(t3,t1) = R2f0(t3,t1) + evops(1,1)%Ufe(j,i,h,f,t3)*extD(i,j,f,h,t1)
                                                               !evops(1,1)%Ufe(i,j,f,h,t3)*extD(i,j,f,h,t1) !W2g(a,b,t3)*D2g(a,b,t1)
                                !end do
                            end do
                        end do
                    end do

                end do
            end do

        else

        end if

    end subroutine make_orR2f


    !
    ! R2f with orientational averaging  SECULAR VERSION
    !
    subroutine make_orR2f_sec(ti,t1)
        integer(i4b) :: ti, t1

        ! local
        integer(i4b) ::  t3, i,j,f,h

        if (ti == 1) then

            !R2f0 = 0.0_dp
            R2f0(:,t1) = 0.0_dp

            !do t1 = 1, Nt1
                do t3 = 1, Nt1

                    do i = 1,Ne
                        do j = 1, Nf
                            !do f = 1, Ne
                                !do h = 1, Nf
                                 !h = 3
                                 h = j
                                 f = i
!                                    R2f0(t3,t1) = R2f0(t3,t1) + evops(1,1)%UfeS(j,i,t3)*extD(i,j,f,h,t1)
                                    R2f0(t3,t1) = R2f0(t3,t1) + evops(1,1)%UfeS(j,i,t3)*extD(i,j,f,h,1)
                                    							! evops(1,1)%Ufe(j,i,h,f,t3)*extD(i,j,f,h,t1)
                                                               !evops(1,1)%Ufe(i,j,f,h,t3)*extD(i,j,f,h,t1) !W2g(a,b,t3)*D2g(a,b,t1)
                                !end do
                            !end do
                        end do
                    end do

                end do
            !end do

        else

        end if

    end subroutine make_orR2f_sec



    !
    ! R3g
    !
    subroutine make_R3g(ti)
        integer :: ti
        ! local
        integer :: t1, t3, a, b

        if (ti == 1) then

            do t1 = 1, Nt1
                do t3 = 1, Nt1
                    R3g0(t3,t1) = 0.0_dp
                    do a = 1, Ng
                        do b = 1, Ng
                            R3g0(t3,t1) = R3g0(t3,t1) + W3g(a,b,t3)*D3g(a,b,t1)
                        end do
                    end do
                end do
            end do

        else

        end if

    end subroutine make_R3g


    !
    ! R3g with orientational averaging
    !
    subroutine make_orR3g(ti)
        integer :: ti

        ! local
        integer :: t1, t3, i,j,h,g

        if (ti == 1) then

            R3g0 = 0.0_dp

            do t1 = 1, Nt1
                do t3 = 1, Nt1

                    do i = 1, Ne
                        do j = 1, Ng
                            do h = 1, Ne
                                do g = 1, Ng
                                    R3g0(t3,t1) = R3g0(t3,t1) + evops(1,1)%Ueg(i,j,h,g,t3)*extD(i,j,h,g,t1) !W2g(a,b,t3)*D2g(a,b,t1)
                                end do
                            end do
                        end do
                    end do

                end do
            end do

        else

        end if

    end subroutine make_orR3g

    !
    ! R3g with orientational averaging SECULAR VERSION
    !
    subroutine make_orR3g_sec(ti)
        integer :: ti

        ! local
        integer :: t1, t3, i,j,h,g

        if (ti == 1) then

            R3g0 = 0.0_dp

            do t1 = 1, Nt1
                do t3 = 1, Nt1

                    do i = 1, Ne
                        do j = 1, Ng
                            !do h = 1, Ne
                                h = i
                                !do g = 1, Ng
                                	g = j
                                    R3g0(t3,t1) = R3g0(t3,t1) + evops(1,1)%Ueg(i,j,h,g,t3)*extD(i,j,h,g,t1) !W2g(a,b,t3)*D2g(a,b,t1)
                                !end do
                            !end do
                        end do
                    end do

                end do
            end do

        else

        end if

    end subroutine make_orR3g_sec


    !
    ! R3f
    !
    subroutine make_R3f(ti)
        integer :: ti
        ! local
        integer :: t1, t3, a, b

        if (ti == 1) then

            do t1 = 1, Nt1
                do t3 = 1, Nt1
                    R3f0(t3,t1) = 0.0_dp
                    do a = 1, Ng
                        do b = 1, Ng
                            R3f0(t3,t1) = R3f0(t3,t1) + W3f(a,b,t3)*D3f(a,b,t1)
                        end do
                    end do
                end do
            end do

        else

        end if

    end subroutine make_R3f


    !
    ! R4g
    !
    subroutine make_R4g(ti)
        integer :: ti
        ! local
        integer :: t1, t3, a, b

        if (ti == 1) then

            do t1 = 1, Nt1
                do t3 = 1, Nt1
                    R4g0(t3,t1) = 0.0_dp
                    do a = 1, Ng
                        do b = 1, Ng
                            R4g0(t3,t1) = R4g0(t3,t1) + W4g(a,b,t3)*D4g(a,b,t1)
                        end do
                    end do
                end do
            end do

        else

        end if

    end subroutine make_R4g


    !
    ! R4g with orientational averaging
    !
    subroutine make_orR4g(ti)
        integer :: ti

        ! local
        integer :: t1, t3, i,j,h,g

        if (ti == 1) then

            R4g0 = 0.0_dp

            do t1 = 1, Nt1
                do t3 = 1, Nt1

                    do i = 1, Ne
                        do j = 1, Ng
                            do h = 1, Ne
                                do g = 1, Ng
                                    R4g0(t3,t1) = R4g0(t3,t1) + evops(1,1)%Ueg(i,j,h,g,t3)*extD(i,j,h,g,t1) !W2g(a,b,t3)*D2g(a,b,t1)
                                end do
                            end do
                        end do
                    end do

                end do
            end do

        else

        end if

    end subroutine make_orR4g

    !
    ! R4g with orientational averaging  SECULAR VERSION
    !
    subroutine make_orR4g_sec(ti)
        integer :: ti

        ! local
        integer :: t1, t3, i,j,h,g

        if (ti == 1) then

            R4g0 = 0.0_dp

            do t1 = 1, Nt1
                do t3 = 1, Nt1

                    do i = 1, Ne
                        do j = 1, Ng
                            !do h = 1, Ne
                            h = i
                                !do g = 1, Ng
                                	g = j
                                    R4g0(t3,t1) = R4g0(t3,t1) + evops(1,1)%Ueg(i,j,h,g,t3)*extD(i,j,h,g,t1) !W2g(a,b,t3)*D2g(a,b,t1)
                                !end do
                            !end do
                        end do
                    end do

                end do
            end do

        else

        end if

    end subroutine make_orR4g_sec


    !
    ! R4f
    !
    subroutine make_R4f(ti)
        integer :: ti
        ! local
        integer :: t1, t3, a, b

        if (ti == 1) then

            do t1 = 1, Nt1
                do t3 = 1, Nt1
                    R4f0(t3,t1) = 0.0_dp
                    do a = 1, Ng
                        do b = 1, Ng
                            R4f0(t3,t1) = R4f0(t3,t1) + W4f(a,b,t3)*D4f(a,b,t1)
                        end do
                    end do
                end do
            end do

        else

        end if

    end subroutine make_R4f


    subroutine test_Upop(U)
        complex(dpc), dimension(:,:,:,:,:), intent(out) :: U

        complex(dpc), dimension(:,:,:,:), allocatable :: unity
        integer :: i,j,k,l,ti

        ALLOCATE(unity,(size(U,1),size(U,2),size(U,3),size(U,4)))


        do i = 1,size(U,1)
            do j = 1, size(U,2)
                do k = 1, size(U,3)
                    do l = 1, size(U,4)
                        if ((i==k).and.(j==l)) then
                            unity(i,j,k,l) = 1.0_dp
                        else
                            unity(i,j,k,l) = 0.0_dp
                        end if
                    end do
                end do
            end do
        end do

        do ti = 1, size(U,5)
            U(:,:,:,:,ti) = unity
        end do

		DEALLOCATE(unity)

    end subroutine test_Upop

    subroutine test_Ucoh(U,oeg,G)
        complex(dpc), dimension(:,:,:,:,:), intent(out) :: U
        real(dp), dimension(:,:), intent(in) :: oeg
        real(dp), intent(in) :: G
        ! local
        !complex(dpc), dimension(:,:,:,:), allocatable :: unity
        integer :: i,j,k,l,ti

        !allocate(unity(size(U,1),size(U,2),size(U,3),size(U,4)))
        do i = 1,size(U,1)
            do j = 1, size(U,2)
                do k = 1, size(U,3)
                    do l = 1, size(U,4)
                        do ti = 1, size(U,5)
                        if ((i==k).and.(j==l)) then
                            U(i,j,k,l,ti) = exp(-(G + (0.0_dp,1.0_dp)*oeg(i,j))*(ti-1)*ddt1)
                        else
                            U(i,j,k,l,ti) = 0.0_dp
                        end if
                        end do
                    end do
                end do
            end do
        end do

    end subroutine test_Ucoh



end module
