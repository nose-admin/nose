module resp3_wd

	use std_types
	use resp3_rss

contains

    !
    ! Window propagators
    !

    !
    ! W1g propagator
    !
    subroutine make_W1g(Win)
        complex(dpc), dimension(:,:,:), intent(out) :: Win

        ! local
        integer :: f,g, h,i,j, ti

        do f = 1, Ne
            do g = 1, Ne
                Win(f,g,:) = 0.0_dp
                do h = 1, Ng
                    do i = 1, Ne
                        do j = 1, Ng
                            do ti = 1, Nt1                 ! dge(j,i)
                                Win(f,g,ti) = Win(f,g,ti) + deg(i,j)*Uegeg(i,j,f,h,ti)*deg(g,h)
                            end do
                        end do
                    end do
                end do
            end do
        end do

    end subroutine make_W1g

    !
    ! W1f propagator
    !
    subroutine make_W1f(Win)
        complex(dpc), dimension(:,:,:), intent(out) :: Win

        ! local
        integer :: f,g, h,i,j, ti

        print *, size(Ufefe,1), size(Ufefe,2), size(Ufefe,3), size(Ufefe,4)

        do f = 1, Ne
            do g = 1, Ne
                Win(f,g,:) = 0.0_dp
                do h = 1, Nf
                    do i = 1, Ne
                        do j = 1, Nf
                            do ti = 1, Nt1                                             ! def(g,h)
                                Win(f,g,ti) = Win(f,g,ti) + dfe(j,i)*conjg(Ufefe(j,i,h,f,ti))*dfe(h,g) !*Uefef(i,j,f,h,ti)*def(g,h)
!                                Win(f,g,ti) = Win(f,g,ti) + dfe(j,i)*conjg(Ufefe(i,j,f,h,ti))*dfe(h,g) !*Uefef(i,j,f,h,ti)*def(g,h)
                            end do
                        end do
                    end do
                end do
            end do
        end do


    end subroutine make_W1f

    !
    ! W3g propagator
    !
    subroutine make_W3g(Win)
        complex(dpc), dimension(:,:,:), intent(out) :: Win

        ! local
        integer :: f,g, h,i,j, ti

        do f = 1, Ng
            do g = 1, Ng
                Win(f,g,:) = 0.0_dp
                do h = 1, Ne
                    do i = 1, Ne
                        do j = 1, Ng
                            do ti = 1, Nt1                 ! dge(j,i)
                                Win(f,g,ti) = Win(f,g,ti) + deg(i,j)*Uegeg(i,j,h,g,ti)*deg(h,f)
                            end do
                        end do
                    end do
                end do
            end do
        end do

    end subroutine make_W3g

    !
    ! W3f propagator
    !
    subroutine make_W3f(Win)
        complex(dpc), dimension(:,:,:), intent(out) :: Win

        ! local
        integer :: f,g, h,i,j, ti

        do f = 1, Ng
            do g = 1, Nf
                Win(f,g,:) = 0.0_dp
                do h = 1, Ne
                    do i = 1, Ne
                        do j = 1, Nf
                            do ti = 1, Nt1                 ! dfe(j,i)
                                Win(f,g,ti) = Win(f,g,ti) + dfe(j,i)*conjg(Ufefe(j,i,g,h,ti))*deg(h,f) !*Uefef(i,j,h,g,ti)*deg(h,f)
                            end do
                        end do
                    end do
                end do
            end do
        end do


    end subroutine make_W3f

    !
    ! W4f propagator
    !
    subroutine make_W4f(Win)
        complex(dpc), dimension(:,:,:), intent(out) :: Win

        ! local
        integer :: f,g, h,i,j, ti

        do f = 1, Nf
            do g = 1, Ng
                Win(f,g,:) = 0.0_dp
                do h = 1, Ne
                    do i = 1, Ne
                        do j = 1, Ng
                            do ti = 1, Nt1                 ! dge(j,i)                   ! def(h,f)
                                Win(f,g,ti) = Win(f,g,ti) + deg(i,j)*Uegeg(i,j,h,g,ti)*dfe(f,h)
                            end do
                        end do
                    end do
                end do
            end do
        end do

    end subroutine make_W4f

    !
    ! Doorway progators
    !

    !
    ! D1g propagator
    !
    subroutine make_D1g(DD)
        complex(dpc), dimension(:,:,:), intent(out) :: DD
        ! local
        integer :: c,e,d,b,a,ti
        do c = 1, Ne
            do e = 1, Ne
                DD(c,e,:) = 0.0_dp
                do d = 1, Ng
                    do b = 1, Ne
                        do a = 1, Ng
                            do ti = 1, Nt1                                         ! dge(d,e)
                                DD(c,e,ti) = DD(c,e,ti) + Uegeg(c,d,b,a,ti)*rhoL(b,a)*deg(e,d)
                            end do
                        end do
                    end do
                end do
            end do
        end do


    end subroutine make_D1g

    !
    ! D1f propagator
    !
    subroutine make_D1f(DD)
        complex(dpc), dimension(:,:,:), intent(out) :: DD

        call make_D1g(DD)

    end subroutine make_D1f

    !
    ! D2g propagator
    !
    subroutine make_D2g(DD)
        complex(dpc), dimension(:,:,:), intent(out) :: DD
        ! local
        integer :: c,e,d,b,a,ti
        do d = 1, Ne
            do e = 1, Ne
                DD(e,d,:) = 0.0_dp
                do c = 1, Ng
                    do b = 1, Ne
                        do a = 1, Ng
                            do ti = 1, Nt1
                                DD(e,d,ti) = DD(e,d,ti) + deg(e,c)*conjg(Uegeg(d,c,b,a,ti))*rhoR(a,b) !Ugege(c,d,a,b,ti)*rhoR(a,b)
                            end do
                        end do
                    end do
                end do
            end do
        end do

    end subroutine make_D2g

    !
    ! D2f propagator
    !
    subroutine make_D2f(DD)
        complex(dpc), dimension(:,:,:), intent(out) :: DD

        call make_D2g(DD)

    end subroutine make_D2f

    !
    ! D3g propagator
    !
    subroutine make_D3g(DD)
        complex(dpc), dimension(:,:,:), intent(out) :: DD
        ! local
        integer :: c,e,d,b,a,ti

        do c = 1, Ng
            do e = 1, Ng
                DD(c,e,:) = 0.0_dp
                do d = 1, Ne
                    do b = 1, Ne
                        do a = 1, Ng
                            do ti = 1, Nt1
                                DD(c,e,ti) = DD(c,e,ti) + conjg(Uegeg(d,c,b,a,ti))*rhoR(a,b)*deg(d,e) !Ugege(c,d,a,b,ti)*rhoR(a,b)*deg(d,e)
                            end do
                        end do
                    end do
                end do
            end do
        end do

    end subroutine make_D3g


    !
    ! D3f propagator
    !
    subroutine make_D3f(DD)
        complex(dpc), dimension(:,:,:), intent(out) :: DD
        ! local
        integer :: c,e,d,b,a,ti

        do c = 1, Ng
            do e = 1, Nf
                DD(c,e,:) = 0.0_dp
                do d = 1, Ne
                    do b = 1, Ne
                        do a = 1, Ng
                            do ti = 1, Nt1
                                DD(c,e,ti) = DD(c,e,ti) + conjg(Uegeg(d,c,b,a,ti))*rhoR(a,b)*dfe(e,d) !Ugege(c,d,a,b,ti)*rhoR(a,b)*def(d,e)
                            end do
                        end do
                    end do
                end do
            end do
        end do

    end subroutine make_D3f

    !
    ! D4g propagator
    !
    subroutine make_D4g(DD)
        complex(dpc), dimension(:,:,:), intent(out) :: DD
        ! local
        integer :: c,e,d,b,a,ti
        do d = 1, Ng
            do e = 1, Ng
                DD(e,d,:) = 0.0_dp
                do c = 1, Ne
                    do b = 1, Ne
                        do a = 1, Ng
                            do ti = 1, Nt1
                                DD(e,d,ti) = DD(e,d,ti) + deg(c,e)*conjg(Uegeg(c,d,b,a,ti))*rhoL(b,a) !Ugege(c,d,a,b,ti)*rhoR(a,b)
                            end do
                        end do
                    end do
                end do
            end do
        end do

    end subroutine make_D4g

    !
    ! D4f propagator
    !
    subroutine make_D4f(DD)
        complex(dpc), dimension(:,:,:), intent(out) :: DD
        ! local
        integer :: c,e,d,b,a,ti
        do d = 1, Ng
            do e = 1, Nf
                DD(e,d,:) = 0.0_dp
                do c = 1, Ne
                    do b = 1, Ne
                        do a = 1, Ng
                            do ti = 1, Nt1                ! def(c,e)
                                DD(e,d,ti) = DD(e,d,ti) + dfe(e,c)*conjg(Uegeg(c,d,b,a,ti))*rhoL(b,a) !Ugege(c,d,a,b,ti)*rhoR(a,b)
                            end do
                        end do
                    end do
                end do
            end do
        end do

    end subroutine make_D4f


end module resp3_wd
