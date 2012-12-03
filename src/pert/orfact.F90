!
!
! Author:  Elizabeth L. Read
!
! Last change: 2007-02-6
!
!
module orfact

    use resources

    implicit none

    real(dp)                                 :: lt(3)
    real(dp), dimension(:,:), allocatable :: GCT
    real(dp), dimension(:,:), allocatable :: GCTS
!   real(dp), dimension(:), allocatable   :: PolA

    integer(dp), dimension(:,:), allocatable :: A2ab_conversion
    integer(dp), dimension(:,:), allocatable :: ab2A_conversion

contains

    !
    ! Initialization
    !
    subroutine init_orfact()

        integer(i4b) :: basis                     ! can be 1 - site or 2 - excitonic
        integer(i4b) :: i, a,b,j,c,k,d,s,p,k1
        real(dp) :: PolV(3,4),trPolV(4,3),LCM(4,4)
        integer(i4b) :: bk

		call print_log_message("Init orfact ...", 8)

        if (allocated(GCT)) deallocate(GCT)
        if (allocated(GCTS)) deallocate(GCTS)

        allocate(GCT(nmax(nr_blocks)+N_ExOrder(nr_blocks,1)* &
                   (N_ExOrder(nr_blocks,2)+1), &
            nmax(nr_blocks)+N_ExOrder(nr_blocks,1)*(N_ExOrder(nr_blocks,2)+1)))
        allocate(GCTS(nmax(nr_blocks)+N_ExOrder(nr_blocks,1)* &
                   (N_ExOrder(nr_blocks,2)+1), &
            nmax(nr_blocks)+N_ExOrder(nr_blocks,1)*(N_ExOrder(nr_blocks,2)+1)))

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
            allocate(ab2A_conversion(k*(k-1)/2,k))
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
            ! excitonic dipole-moments
            !

            ! One-excitons
            !
            do i = 1, N_ExOrder(bk,1)
                mun(i,1) = current_e_block%dx(i)
                mun(i,2) = current_e_block%dy(i)
                mun(i,3) = current_e_block%dz(i)
                !print *, mun(i,1)
            end do
            !stop

            !
            ! two-excitons
            !
            if (use_twoexcitons) then
            ! transition from an exciton j to a double-exciton (s,p)
            k = 0
            !i = 0
            do i = 1, N_ExOrder(bk,2) !N1
                ! two-exciton i
                do j = 1, N_ExOrder(bk,1)
                    ! one-exciton j

                    k = k + 1

                    mun(N_ExOrder(bk,1)+k,1) = current_e_block%dx_2(i,j)
                    mun(N_ExOrder(bk,1)+k,2) = current_e_block%dy_2(i,j)
                    mun(N_ExOrder(bk,1)+k,3) = current_e_block%dz_2(i,j)

                end do
            !end do
            end do
            end if

            !
            ! site basis dipole-moments
            !

            ! One excitation
            !
            do i = 1, N_ExOrder(current_s_block%id,1)
                munS(i,1) = current_s_block%dx(i)
                munS(i,2) = current_s_block%dy(i)
                munS(i,3) = current_s_block%dz(i)
            end do

            ! Two excitations
            !
            ! transition from an exciton j to a double-exciton (s,p)
            k = 0
            i = 0
            do s = 1, N_ExOrder(bk,1) !N1
            do p = s+1, N_ExOrder(bk,1) !N1
                i = i + 1
                ! two-exciton (s,p)
                do j = 1, N_ExOrder(bk,1)
                    ! one-exciton j

                    k = k + 1

                    if (s == j) then
                        munS(N_ExOrder(bk,1)+k,1) = current_s_block%dx(p)
                        munS(N_ExOrder(bk,1)+k,2) = current_s_block%dy(p)
                        munS(N_ExOrder(bk,1)+k,3) = current_s_block%dz(p)
                    else if (p == j) then
                        munS(N_ExOrder(bk,1)+k,1) = current_s_block%dx(s)
                        munS(N_ExOrder(bk,1)+k,2) = current_s_block%dy(s)
                        munS(N_ExOrder(bk,1)+k,3) = current_s_block%dz(s)
                    else
                        munS(N_ExOrder(bk,1)+k,1) = 0.0_dp
                        munS(N_ExOrder(bk,1)+k,2) = 0.0_dp
                        munS(N_ExOrder(bk,1)+k,3) = 0.0_dp
                    end if
                end do
            end do
            end do



            !
            ! Cosine table for excitons
            !

            !Greek Cosine Table
            do a=1,1+N_ExOrder(bk,2)
            do b=1,N_ExOrder(bk,1)
                if (a == 1) then
                    j = b
                else
                    j = N_ExOrder(bk,1)+ab2A_conversion(a-1,b)
                end if
                do c=1,1+N_ExOrder(bk,2)
                do d=1,N_ExOrder(bk,1)
                    if (c == 1) then
                        k = d
                    else
                        k = N_ExOrder(bk,1)+ab2A_conversion(c-1,d)
                    end if !d !nmax(i) + (c-1)*N_ExOrder(i,1) + d

                    GCT(j,k)=dot_product(mun(j,:),mun(k,:))

                end do
                end do
            end do
            end do

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

            !
            ! Cosine table for site basis
            !
            do a=1,1+N_ExOrder(bk,2)
            do b=1,N_ExOrder(bk,1)
                if (a == 1) then
                    j = b
                else
                    j = N_ExOrder(bk,1)+ab2A_conversion(a-1,b)
                end if
                do c=1,1+N_ExOrder(bk,2)
                do d=1,N_ExOrder(bk,1)
                    if (c == 1) then
                        k = d
                    else
                        k = N_ExOrder(bk,1)+ab2A_conversion(c-1,d)
                    end if !d !nmax(i) + (c-1)*N_ExOrder(i,1) + d

                    GCTS(j,k)=dot_product(munS(j,:),munS(k,:))

                end do
                end do
            end do
            end do



            if (.not.resources_have_next_block()) exit
            call resources_next_block()
            bk = bk + 1


        end do

		call print_log_message("... orfact done.", 8)


    end subroutine init_orfact

    !
    ! Cleaning the module
    !
    subroutine clean_orfact()

        if (allocated(GCT)) then
            deallocate(GCT)
        end if
        if (allocated(GCTS)) then
            deallocate(GCTS)
        end if

    end subroutine clean_orfact

    !
    ! Orientational factors of Liouville pathways
    !

    !
    !
    !R2,R3
    function orfact_23(ea,na,eb,nb,ec,nc,ed,nd) result(ORF)
        integer(i4b), intent(in) :: ea,eb,ec,ed,na,nb,nc,nd
        integer(i4b) :: ga,gb,gc,gd
        real(dp) :: ORF,gtt(3)

        ga = nmax(na)+ea
        gb = nmax(nb)+eb
        gc = nmax(nc)+ec
        gd = nmax(nd)+ed

        gtt(1) = GCT(ga,gb)*GCT(gc,gd)
        gtt(2) = GCT(ga,gc)*GCT(gb,gd)
        gtt(3) = GCT(ga,gd)*GCT(gb,gc)

        ORF = (1.0_dp/30.0_dp)*(gtt(1)*lt(1)+gtt(2)*lt(2)+gtt(3)*lt(3))

    end function orfact_23

    !
    ! Pathways including ground and one-exciton states / site basis
    !
    function orfact_23_S(ea,na,eb,nb,ec,nc,ed,nd) result(ORF)
        integer(i4b), intent(in) :: ea,eb,ec,ed,na,nb,nc,nd
        integer(i4b) :: ga,gb,gc,gd
        real(dp) :: ORF,gtt(3)

        logical :: switch

        switch = use_twoexcitons  ! here we should have a more specific switch (use_excitons or so)

        if (switch) then

            ORF = orfact_23(ea,na,eb,nb,ec,nc,ed,nd)

        else

        ga = nmax(na)+ea
        gb = nmax(nb)+eb
        gc = nmax(nc)+ec
        gd = nmax(nd)+ed

        gtt(1) = GCTS(ga,gb)*GCTS(gc,gd)
        gtt(2) = GCTS(ga,gc)*GCTS(gb,gd)
        gtt(3) = GCTS(ga,gd)*GCTS(gb,gc)

        ORF = (1.0_dp/30.0_dp)*(gtt(1)*lt(1)+gtt(2)*lt(2)+gtt(3)*lt(3))

        end if

    end function orfact_23_S



    !R1
    function orfact_1(f,nf,ea,eb,ec,nc,ed,nd) result(ORF)
        integer(i4b) :: ga,gb,gc,gd
        integer(i4b), intent(in) :: f,ea,eb,ec,ed,nf,nc,nd
        real(dp) :: ORF,gtt(3)

        ga = nmax(nf)+f*N_ExOrder(nf,1)+ea
        gb = nmax(nf)+f*N_ExOrder(nf,1)+eb
        gc = nmax(nc)+ec
        gd = nmax(nd)+ed

        gtt(1) = GCT(ga,gb)*GCT(gc,gd)
        gtt(2) = GCT(ga,gc)*GCT(gb,gd)
        gtt(3) = GCT(ga,gd)*GCT(gb,gc)

        ORF = (1.0_dp/30.0_dp)*(gtt(1)*lt(1)+gtt(2)*lt(2)+gtt(3)*lt(3))

    end function orfact_1

    !
    ! Pathways including doubly-excited states / excitonic basis
    !
    function orfact_e(g,ng,ea,na,f,nf,eb,nb,ec,nc,ed,nd) result(ORF)
        integer(i4b) :: ga,gb,gc,gd,fa,fb
        integer(i4b), intent(in) :: g,ng,ea,na,f,nf,eb,nb,ec,nc,ed,nd
        real(dp) :: ORF,gtt(3)

        fa = N_ExOrder(ng,1) + ab2A_conversion(g,ea)
        fb = N_ExOrder(nf,1) + ab2A_conversion(f,eb)
        ga = nmax(ng)+fa !g*N_ExOrder(ng,1)+ea
        gb = nmax(nf)+fb !f*N_ExOrder(nf,1)+eb
        gc = nmax(nc)+ec
        gd = nmax(nd)+ed

        gtt(1) = GCT(ga,gb)*GCT(gc,gd)
        gtt(2) = GCT(ga,gc)*GCT(gb,gd)
        gtt(3) = GCT(ga,gd)*GCT(gb,gc)

        ORF = (1.0_dp/30.0_dp)*(gtt(1)*lt(1)+gtt(2)*lt(2)+gtt(3)*lt(3))


    end function orfact_e


    !
    ! Pathways including doubly-excited states / site basis
    !
    function orfact_e_S(g,ng,ea,na,f,nf,eb,nb,ec,nc,ed,nd) result(ORF)
        integer(i4b) :: ga,gb,gc,gd
        integer(i4b), intent(in) :: g,ng,ea,na,f,nf,eb,nb,ec,nc,ed,nd
        real(dp) :: ORF,gtt(3)
        integer(i4b) :: fa, fb

        logical :: switch

        switch = use_twoexcitons  ! here we should have a more specific switch (use_excitons or so)

        if (switch) then

            ORF = orfact_e(g,ng,ea,na,f,nf,eb,nb,ec,nc,ed,nd)

        else

        fa = N_ExOrder(ng,1) + ab2A_conversion(g,ea)
        fb = N_ExOrder(nf,1) + ab2A_conversion(f,eb)
        ga = nmax(ng) + fa !g*N_ExOrder(ng,1)+ea
        gb = nmax(nf) + fb !f*N_ExOrder(nf,1)+eb
        gc = nmax(nc) + ec
        gd = nmax(nd) + ed

        gtt(1) = GCTS(ga,gb)*GCTS(gc,gd)
        gtt(2) = GCTS(ga,gc)*GCTS(gb,gd)
        gtt(3) = GCTS(ga,gd)*GCTS(gb,gc)

        ORF = (1.0_dp/30.0_dp)*(gtt(1)*lt(1)+gtt(2)*lt(2)+gtt(3)*lt(3))

        end if

    end function orfact_e_S


end module orfact
