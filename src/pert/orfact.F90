#include "util_allocation.h"
!
!
! Author:  Elizabeth L. Read
!
! Last change: 2012-04-24 (Jan Olsina)
!
!
module orfact

    use resources

    implicit none

    real(dp)                                 :: lt(3)
    real(dp), dimension(:,:,:,:), allocatable, private :: GCTegeg
    real(dp), dimension(:,:,:,:), allocatable, private :: GCTSegeg
    real(dp), dimension(:,:,:,:), allocatable, private :: GCTfeeg
    real(dp), dimension(:,:,:,:), allocatable, private :: GCTSfeeg
    real(dp), dimension(:,:,:,:), allocatable, private :: GCTfefe
    real(dp), dimension(:,:,:,:), allocatable, private :: GCTSfefe
!   real(dp), dimension(:), allocatable   :: PolA

    integer(dp), dimension(:,:), allocatable, private :: A2ab_conversion
    integer(dp), dimension(:,:), allocatable, private :: ab2A_conversion

    logical, private :: orfact_initiated = .false.
    integer(i4b), private :: Ng, Ne, Nf

contains

    !
    ! Initialization
    !
    subroutine init_orfact()

        integer(i4b) :: basis                     ! can be 1 - site or 2 - excitonic
        integer(i4b) :: i, a,b,j,c,k,d,s,p,k1,l
        real(dp) :: PolV(3,4),trPolV(4,3),LCM(4,4)
        integer(i4b) :: bk
        real(dp), dimension(:,:), allocatable :: dx2,dy2,dz2

		call print_log_message("Init orfact ...", 8)
		write(*,*) 'orfact',allocated(GCTegeg),allocated(GCTSegeg),allocated(GCTfeeg),allocated(GCTSfeeg),allocated(GCTfefe),allocated(GCTSfefe)
		write(*,*) 'orfact',associated(evops)
		if(associated(evops)) then
			write(*,*) 'orfact',associated(evops(1,1)%Ueg),associated(evops(1,1)%Uee),associated(evops(1,1)%Ufe)
		endif
		call flush()

		if(.not. associated(evops)) then
			return
			!call print_error_message(-1,"orfact initiated before evops are allocated!")
		end if

        if (allocated(GCTegeg)) then
        	DEALLOCATE(GCTegeg)
        end if
        if (allocated(GCTSegeg)) then
    	    DEALLOCATE(GCTSegeg)
        end if
        if (allocated(GCTfeeg)) then
	        DEALLOCATE(GCTfeeg)
        end if
        if (allocated(GCTSfeeg)) then
        	DEALLOCATE(GCTSfeeg)
        end if
        if (allocated(GCTfefe)) then
    	    DEALLOCATE(GCTfefe)
        end if
        if (allocated(GCTSfefe)) then
	        DEALLOCATE(GCTSfefe)
        end if

        Ng = size(evops(1,1)%Ueg,2)
        Ne = size(evops(1,1)%Ueg,1)
        if(associated(evops(1,1)%Ufe)) then
        	Nf = size(evops(1,1)%Ufe,1)
        elseif(associated(evops(1,1)%UfeS) ) then
        	Nf = size(evops(1,1)%UfeS,1)
        else
        	Nf = 0
        end if

        if(.not. use_twoexcitons) then
        	Nf = 0
        end if

        write(*,*) 'Nf=',Nf

        ALLOCATE(GCTegeg,(Ne,Ng,Ne,Ng))
        ALLOCATE(GCTSegeg,(Ne,Ng,Ne,Ng))

        if(Nf /= 0) then
        	ALLOCATE(GCTfeeg,(Nf,Ne,Ne,Ng))
        	ALLOCATE(GCTSfeeg,(Nf,Ne,Ne,Ng))
        	ALLOCATE(GCTfefe,(Nf,Ne,Nf,Ne))
        	ALLOCATE(GCTSfefe,(Nf,Ne,Nf,Ne))
        end if

        do i=1,4
            PolV(1,i) = 0.0_dp
            PolV(2,i) = SIN(PolA(i))
            PolV(3,i) = COS(PolA(i))
        end do

        ! Calculate three latin index factors determined by field
        trPolV=transpose(PolV)

        !Latin Cosine Table
        LCM=matmul(trPolV,PolV)
        lt(1)=4.0_dp*LCM(1,2)*LCM(3,4)-LCM(1,3)*LCM(2,4)-LCM(1,4)*LCM(2,3)
        lt(2)=4.0_dp*LCM(1,3)*LCM(2,4)-LCM(1,2)*LCM(3,4)-LCM(1,4)*LCM(2,3)
        lt(3)=4.0_dp*LCM(1,4)*LCM(2,3)-LCM(1,2)*LCM(3,4)-LCM(1,3)*LCM(2,4)

        k = maxval(Block_Sizes)
        if (.not.allocated(ab2A_conversion)) then
            ALLOCATE(ab2A_conversion,(k*(k-1)/2,k))
        end if



        call resources_rewind_blocks()
        bk = 1
        do


            ! index conversion in a block
            k = 0
            do i = 1, N_ExOrder(bk,2) !N1
                ! two-exciton i
                do j = 1, N_ExOrder(bk,1)
                    ! one-exciton j
                    k = k + 1
                    ab2A_conversion(i,j) = k
                end do
            end do
!
!
!            !
!            ! excitonic dipole-moments
!            !
!
!            ! One-excitons
!            !
!            do i = 1, N_ExOrder(bk,1)
!                mun(i,1) = current_e_block%dx(i)
!                mun(i,2) = current_e_block%dy(i)
!                mun(i,3) = current_e_block%dz(i)
!                !print *, mun(i,1)
!            end do
!            !stop
!
!            !
!            ! two-excitons
!            !
!            if (use_twoexcitons) then
!            ! transition from an exciton j to a double-exciton (s,p)
!            k = 0
!            !i = 0
!            do i = 1, N_ExOrder(bk,2) !N1
!                ! two-exciton i
!                do j = 1, N_ExOrder(bk,1)
!                    ! one-exciton j
!
!                    k = k + 1
!
!                    mun(N_ExOrder(bk,1)+k,1) = current_e_block%dx_2(i,j)
!                    mun(N_ExOrder(bk,1)+k,2) = current_e_block%dy_2(i,j)
!                    mun(N_ExOrder(bk,1)+k,3) = current_e_block%dz_2(i,j)
!
!                end do
!            !end do
!            end do
!            end if
!
!            !
!            ! site basis dipole-moments
!            !
!
!            ! One excitation
!            !
!            do i = 1, N_ExOrder(current_s_block%id,1)
!                munS(i,1) = current_s_block%dx(i)
!                munS(i,2) = current_s_block%dy(i)
!                munS(i,3) = current_s_block%dz(i)
!            end do
!
!            ! Two excitations
!            !
!            ! transition from an exciton j to a double-exciton (s,p)
!            k = 0
!            i = 0
!            do s = 1, N_ExOrder(bk,1) !N1
!            do p = s+1, N_ExOrder(bk,1) !N1
!                i = i + 1
!                ! two-exciton (s,p)
!                do j = 1, N_ExOrder(bk,1)
!                    ! one-exciton j
!
!                    k = k + 1
!
!                    if (s == j) then
!                        munS(N_ExOrder(bk,1)+k,1) = current_s_block%dx(p)
!                        munS(N_ExOrder(bk,1)+k,2) = current_s_block%dy(p)
!                        munS(N_ExOrder(bk,1)+k,3) = current_s_block%dz(p)
!                    else if (p == j) then
!                        munS(N_ExOrder(bk,1)+k,1) = current_s_block%dx(s)
!                        munS(N_ExOrder(bk,1)+k,2) = current_s_block%dy(s)
!                        munS(N_ExOrder(bk,1)+k,3) = current_s_block%dz(s)
!                    else
!                        munS(N_ExOrder(bk,1)+k,1) = 0.0_dp
!                        munS(N_ExOrder(bk,1)+k,2) = 0.0_dp
!                        munS(N_ExOrder(bk,1)+k,3) = 0.0_dp
!                    end if
!                end do
!            end do
!            end do
!
!
!
!            !
!            ! Cosine table for excitons
!            !
!
!            !Greek Cosine Table
!            do a=1,1+N_ExOrder(bk,2)
!            do b=1,N_ExOrder(bk,1)
!                if (a == 1) then
!                    j = b
!                else
!                    j = N_ExOrder(bk,1)+ab2A_conversion(a-1,b)
!                end if
!                do c=1,1+N_ExOrder(bk,2)
!                do d=1,N_ExOrder(bk,1)
!                    if (c == 1) then
!                        k = d
!                    else
!                        k = N_ExOrder(bk,1)+ab2A_conversion(c-1,d)
!                    end if !d !nmax(i) + (c-1)*N_ExOrder(i,1) + d
!
!                    GCT(j,k)=dot_product(mun(j,:),mun(k,:))
!
!                end do
!                end do
!            end do
!            end do

!        do a=1,N_ExOrder(i,2)+1
!        do b=1,N_ExOrder(i,1)
!            j = nmax(i) + (a-1)*N_ExOrder(i,1) + b
!            do c=1,N_ExOrder(i,2)+1
!            do d=1,N_ExOrder(i,1)
!               k = nmax(i) + (c-1)*N_ExOrder(i,1) + d
!               !k1 = nmax(i)+ N_ExOrder(i,1) + ab2A_conversion(g,ea)
!               GCT(j,k)=dot_product(mun(j,:),mun(k,:))
!               !print *, j,k,k1, GCT(j,k)
!            end do
!            end do
!        end do
!        end do

!        do a=1,N_ExOrder(i,2)+1
!        do b=1,N_ExOrder(i,1)
!            j = nmax(i) + (a-1)*N_ExOrder(i,1) + b
!            do c=1,N_ExOrder(i,2)+1
!            do d=1,N_ExOrder(i,1)
!               k = nmax(i) + (c-1)*N_ExOrder(i,1) + d
!               !k1 = nmax(i)+ N_ExOrder(i,1) + ab2A_conversion(g,ea)
!               GCTS(j,k)=dot_product(munS(j,:),munS(k,:))
!               !print *, j,k,k1, GCT(j,k)
!            end do
!            end do
!        end do
!        end do

!            !
!            ! Cosine table for site basis
!            !
!            do a=1,1+N_ExOrder(bk,2)
!            do b=1,N_ExOrder(bk,1)
!                if (a == 1) then
!                    j = b
!                else
!                    j = N_ExOrder(bk,1)+ab2A_conversion(a-1,b)
!                end if
!                do c=1,1+N_ExOrder(bk,2)
!                do d=1,N_ExOrder(bk,1)
!                    if (c == 1) then
!                        k = d
!                    else
!                        k = N_ExOrder(bk,1)+ab2A_conversion(c-1,d)
!                    end if !d !nmax(i) + (c-1)*N_ExOrder(i,1) + d
!
!                    GCTS(j,k)=dot_product(munS(j,:),munS(k,:))
!
!                end do
!                end do
!            end do
!            end do

            if (use_twoexcitons .and. (.not. special_carotenoid_hack)) then
                ALLOCATE(dx2,(Nf,Ne))
                ALLOCATE(dy2,(Nf,Ne))
                ALLOCATE(dz2,(Nf,Ne))

                do i = 1, Nf
                    k = current_e_block%ione1(i)
                    l = current_e_block%ione2(i)
                    do j = 1, Ne
                        if (j == k) then
                            dx2(i,j) = current_s_block%dx(l,1)
                            dy2(i,j) = current_s_block%dy(l,1)
                            dz2(i,j) = current_s_block%dz(l,1)
                        else if (j == l) then
                            dx2(i,j) = current_s_block%dx(k,1)
                            dy2(i,j) = current_s_block%dy(k,1)
                            dz2(i,j) = current_s_block%dz(k,1)
                        else
                            dx2(i,j) = 0.0_dp
                            dy2(i,j) = 0.0_dp
                            dz2(i,j) = 0.0_dp
                        end if
                    end do
                end do
            end if

			!
			!  Greek Cosine Tables
			!
			do i=1,Ne
			do j=1,Ng
			do k=1,Ne
			do l=1,Ng
				GCTegeg(i,j,k,l) =										&
					current_e_block%dx(i,j)*current_e_block%dx(k,l) + 	&
					current_e_block%dy(i,j)*current_e_block%dy(k,l) + 	&
					current_e_block%dz(i,j)*current_e_block%dz(k,l)

				GCTSegeg(i,j,k,l) =										&
					current_s_block%dx(i,j)*current_s_block%dx(k,l) + 	&
					current_s_block%dy(i,j)*current_s_block%dy(k,l) + 	&
					current_s_block%dz(i,j)*current_s_block%dz(k,l)
			end do
			end do
			end do
			end do

			if(Nf /= 0 .and. use_twoexcitons) then

			do i=1,Nf
			do j=1,Ne
			do k=1,Ne
			do l=1,Ng
				GCTfeeg(i,j,k,l) = 											&
					current_e_block%dx_2(i,j)*current_e_block%dx(k,l) + 	&
					current_e_block%dy_2(i,j)*current_e_block%dy(k,l) + 	&
					current_e_block%dz_2(i,j)*current_e_block%dz(k,l)

			if(special_carotenoid_hack) then
				GCTSfeeg(i,j,k,l) = GCTfeeg(i,j,k,l)
			else
				GCTSfeeg(i,j,k,l) = 										&
					dx2(i,j)*current_s_block%dx(k,l) + 	&
					dy2(i,j)*current_s_block%dy(k,l) + 	&
					dz2(i,j)*current_s_block%dz(k,l)
			end if

			end do
			end do
			end do
			end do

			do i=1,Nf
			do j=1,Ne
			do k=1,Nf
			do l=1,Ne
				GCTfefe(i,j,k,l) = 											&
					current_e_block%dx_2(i,j)*current_e_block%dx_2(k,l) + 	&
					current_e_block%dy_2(i,j)*current_e_block%dy_2(k,l) + 	&
					current_e_block%dz_2(i,j)*current_e_block%dz_2(k,l)

			if(special_carotenoid_hack) then
				GCTSfefe(i,j,k,l) = GCTfefe(i,j,k,l)
			else
				GCTSfefe(i,j,k,l) = 										&
					dx2(i,j)*dx2(k,l) + 	&
					dy2(i,j)*dy2(k,l) + 	&
					dz2(i,j)*dz2(k,l)
			end if

			end do
			end do
			end do
			end do

			end if



            if (.not.resources_have_next_block()) exit
            call resources_next_block()
            bk = bk + 1


        end do

		call print_log_message("... orfact done.", 8)

		orfact_initiated = .true.


    end subroutine init_orfact

    !
    ! Cleaning the module
    !
    subroutine clean_orfact()

        if (allocated(GCTegeg)) then
            DEALLOCATE(GCTegeg)
        end if
        if (allocated(GCTSegeg)) then
            DEALLOCATE(GCTSegeg)
        end if
        if (allocated(GCTfeeg)) then
            DEALLOCATE(GCTfeeg)
        end if
        if (allocated(GCTSfeeg)) then
            DEALLOCATE(GCTSfeeg)
        end if
        if (allocated(GCTfefe)) then
            DEALLOCATE(GCTfefe)
        end if
        if (allocated(GCTSfefe)) then
            DEALLOCATE(GCTSfefe)
        end if

        orfact_initiated = .false.

    end subroutine clean_orfact

    !
    ! Orientational factors of Liouville pathways
    !

    !
    !
    !R2,R3
    function orfact_23(ea,ga,eb,gb,ec,gc,ed,gd) result(ORF)
        integer(i4b), intent(in) :: ea,ga,eb,gb,ec,gc,ed,gd
!        integer(i4b) :: gga,ggb,ggc,ggd,eea,eeb,eec,eed
        real(dp) :: ORF,gtt(3)

        gtt(1) = GCTegeg(ea,ga,eb,gb)*GCTegeg(ec,gc,ed,gd)
        gtt(2) = GCTegeg(ea,ga,ec,gc)*GCTegeg(eb,gb,ed,gd)
        gtt(3) = GCTegeg(ea,ga,ed,gd)*GCTegeg(eb,gb,ec,gc)

        ORF = (1.0_dp/30.0_dp)*(gtt(1)*lt(1)+gtt(2)*lt(2)+gtt(3)*lt(3))

    end function orfact_23

    !
    ! Pathways including ground and one-exciton states / site basis
    !
    function orfact_23_S(ea,ga,eb,gb,ec,gc,ed,gd) result(ORF)
        integer(i4b), intent(in) :: ea,ga,eb,gb,ec,gc,ed,gd
!        integer(i4b) :: gga,ggb,ggc,ggd,eea,eeb,eec,eed
        real(dp) :: ORF,gtt(3)

        logical :: switch

        switch = use_twoexcitons  ! here we should have a more specific switch (use_excitons or so)

        if (switch) then

            ORF = orfact_23(ea,ga,eb,gb,ec,gc,ed,gd)

        else

!        eea = nmax(na)+ea
!        eeb = nmax(nb)+eb
!        eec = nmax(nc)+ec
!        eed = nmax(nd)+ed
!        gga = nmax(na)+ga
!        ggb = nmax(nb)+gb
!        ggc = nmax(nc)+gc
!        ggd = nmax(nd)+gd

        gtt(1) = GCTSegeg(ea,ga,eb,gb)*GCTSegeg(ec,gc,ed,gd)
        gtt(2) = GCTSegeg(ea,ga,ec,gc)*GCTSegeg(eb,gb,ed,gd)
        gtt(3) = GCTSegeg(ea,ga,ed,gd)*GCTSegeg(eb,gb,ec,gc)

        ORF = (1.0_dp/30.0_dp)*(gtt(1)*lt(1)+gtt(2)*lt(2)+gtt(3)*lt(3))

        end if

    end function orfact_23_S



    !R1
    function orfact_1(f,ea,eb,gb,ec,gc,ed,gd) result(ORF)
!        integer(i4b) :: ga,gb,gc,gd
        integer(i4b), intent(in) :: f,ea,eb,gb,ec,gc,ed,gd
        real(dp) :: ORF,gtt(3)

!        ga = nmax(nf)+f*N_ExOrder(nf,1)+ea
!        gb = nmax(nf)+f*N_ExOrder(nf,1)+eb
!        gc = nmax(nc)+ec
!        gd = nmax(nd)+ed

        gtt(1) = GCTfeeg(f ,ea,     eb,gb)*GCTegeg(ec,gc,     ed,gd)
        gtt(2) = GCTfeeg(f ,ea,     ec,gc)*GCTegeg(eb,gb,     ed,gd)
        gtt(3) = GCTfeeg(f ,ea,     ed,gd)*GCTegeg(eb,gb,     ec,gc)

        ORF = (1.0_dp/30.0_dp)*(gtt(1)*lt(1)+gtt(2)*lt(2)+gtt(3)*lt(3))

    end function orfact_1

    !
    ! Pathways including doubly-excited states / excitonic basis
    !ea,ga,na,eb,gb,nb,ec,gc,nc,ed,gd,nd
    function orfact_e(g,ea,f,eb,ec,gc,ed,gd) result(ORF)
!        integer(i4b) :: ga,gb,gc,gd,fa,fb
        integer(i4b), intent(in) :: g,ea,f,eb,ec,gc,ed,gd
        real(dp) :: ORF,gtt(3)

!        fa = N_ExOrder(ng,1) + ab2A_conversion(g,ea)
!        fb = N_ExOrder(nf,1) + ab2A_conversion(f,eb)
!        ga = nmax(ng)+fa !g*N_ExOrder(ng,1)+ea
!        gb = nmax(nf)+fb !f*N_ExOrder(nf,1)+eb
!        gc = nmax(nc)+ec
!        gd = nmax(nd)+ed

        gtt(1) = GCTfefe(g,ea,    f ,eb)*GCTegeg(ec,gc,    ed,gd)
        gtt(2) = GCTfeeg(g,ea,    ec,gc)*GCTfeeg(f ,eb,    ed,gd)
        gtt(3) = GCTfeeg(g,ea,    ed,gd)*GCTfeeg(f ,eb,    ec,gc)

        ORF = (1.0_dp/30.0_dp)*(gtt(1)*lt(1)+gtt(2)*lt(2)+gtt(3)*lt(3))


    end function orfact_e


    !
    ! Pathways including doubly-excited states / site basis
    !
    function orfact_e_S(g,ea,f,eb,ec,gc,ed,gd) result(ORF)
!        integer(i4b) :: ga,gb,gc,gd
        integer(i4b), intent(in) :: g,ea,f,eb,ec,gc,ed,gd
        real(dp) :: ORF,gtt(3)
        integer(i4b) :: fa, fb

        logical :: switch


        switch = use_twoexcitons  ! here we should have a more specific switch (use_excitons or so)

        if (switch) then

            ORF = orfact_e(g,ea,f,eb,ec,gc,ed,gd)

        else

!        fa = N_ExOrder(ng,1) + ab2A_conversion(g,ea)
!        fb = N_ExOrder(nf,1) + ab2A_conversion(f,eb)
!        ga = nmax(ng) + fa !g*N_ExOrder(ng,1)+ea
!        gb = nmax(nf) + fb !f*N_ExOrder(nf,1)+eb
!        gc = nmax(nc) + ec
!        gd = nmax(nd) + ed

        gtt(1) = GCTfefe(g,ea,    f ,eb)*GCTegeg(ec,gc,    ed,gd)
        gtt(2) = GCTfeeg(g,ea,    ec,gc)*GCTfeeg(f ,eb,    ed,gd)
        gtt(3) = GCTfeeg(g,ea,    ed,gd)*GCTfeeg(f ,eb,    ec,gc)

        ORF = (1.0_dp/30.0_dp)*(gtt(1)*lt(1)+gtt(2)*lt(2)+gtt(3)*lt(3))

        end if

    end function orfact_e_S


end module orfact
