!
!
!
!
module qme_weak_excitons

    use prepare
	use twod

	use resources_qme
	use std_types

	use numer_interp
	use numer_ode
	use numer_fft
	use sci_misc

	implicit none

	complex(dpc), dimension(:,:,:,:,:), allocatable :: R1, R2 ! for the relaxation matrix

contains


!------------------------------------- "create" functions -----------------------------------------


	!
	! Calculating evolution of populations (by a simple rate equation)
	!
	subroutine weak_excitons_pop()

		integer :: kk,k,i
		real(dp) :: t
		real(dp), dimension(:), allocatable :: rho, drho, dd
		real(dp) :: Nrm

		call resources_rewind_blocks()

		kk = 1
		! loop over blocks
		do

			allocate(dd(N1))
			Nrm = 0.0_dp
			do i = 1, N1
			   dd(i) = dx(i,1)*dx(i,1) + dy(i,1)*dy(i,1) + dz(i,1)*dz(i,1)
			   Nrm = Nrm + dd(i)	   ! Normalisation
			end do
			allocate(rho(N1),drho(N1))

			! initial condition set by the transition dipoles
			rho = dd/Nrm

			call create_rates()

			i = 0
			! loop over time
			do
				i = i + 1
				t = (i-1)*gt(1)*dt

				! save the populations
				do k = kk, kk+N1-1
					pops%P(i,k) = rho(k-kk+1)
				end do

			   if (i == Nt(1)) exit

			  call make_step_pop(t,rho,drho)

			end do


			deallocate(rho,drho,dd)

			kk = kk + N1

			if (.not.resources_have_next_block()) exit

			call resources_next_block()

		end do

	end subroutine weak_excitons_pop


	!
	! Calculating evolution of the optical coherences
	!
	subroutine weak_excitons_coh()

		integer :: kk,k,i, j
		real(dp) :: t
		complex(dpc), dimension(:), allocatable :: gc, dgc

		call resources_rewind_blocks()

		kk = 1
		! loop over blocks
		do

			allocate(cpl(Nt(1),N1,N1)) ! prepapring matrix of couplings
			cpl = 0.0_dp
			call create_couplings()
			allocate(RelMa(N1,N1,N1,N1,Nt(1))) ! preparing relaxation matrix
			RelMa = 0.0_dp
			call create_relaxation()

			! now its allocated here, but will be
			evops(kk,kk)%Ueg = 0.0_dp

			allocate(gc(N1),dgc(N1))

			! loop over initial conditions
			do j = 1,N1

			gc = 0.0_dp
			gc(j) = 1.0_dp !dd

			i = 0
			! loop over time
			do

				i = i + 1
				t = (i-1)*gt(1)*dt

				! temporary - no evolution due to coupling
				!gc = 1.0_dp

				! save the coherences !IMPORTANT: switching back from interaction picture coherences here
				do k = kk, kk+N1-1			!down here, from *
					gcohs%C(i,k) = gc(k-kk+1)*exp(-all_goft(current_s_block%gindex(k-kk+1))%gg(gt(1)*(i-1)+1)&
					-(0.0_dp,1.0_dp)*(current_s_block%en(k-kk+1)-rwa)*t)
					evops(kk,kk)%Ueg(k-kk+1,1,j,1,i) = gcohs%C(i,k)
				end do

				if (i == Nt(1)) exit

				call make_step_coh(t,gc,dgc)

			end do

			end do

			deallocate(gc,dgc,cpl,RelMa)

			kk = kk + N1

		   if (.not.resources_have_next_block()) exit

			call resources_next_block()

		end do

	end subroutine weak_excitons_coh


	!
	! Calculating evolution of the non-optical coherences + corresponding populations
	!
	subroutine weak_excitons_pop_coh()

		integer  :: kk,k,j,i,a,b
		real(dp) :: t
		complex(dpc), dimension(:,:), allocatable :: prc, dprc
		real(dp), dimension(:), allocatable :: dd


		call resources_rewind_blocks()

		kk = 1

		! loop over blocks
		do

			allocate(prc(N1,N1),dprc(N1,N1))

			allocate(cpl(Nt(1),N1,N1)) ! preparing matrix of couplings
			cpl = 0.0_dp
			call create_couplings()
			allocate(RelMa(N1,N1,N1,N1,Nt(1))) ! preparing relaxation matrix
			RelMa = 0.0_dp
			call create_relaxation()

			evops(kk,kk)%Uee = 0.0_dp


		! loop over i.c.'s
		do b = 1, N1
		  do a = 1, N1

			  prc = 0.0_dp; dprc = 0.0_dp
			  prc(a,b) = 1.0_dp

			i = 0
			! loop over time
			do

				i = i + 1
				t = (i-1)*gt(1)*dt

				! save the (N1 x N1) matrix elements into the (time x N1 x N1) matrix of time evolution
				do j = kk, kk+N1-1	  ! because row index (k here) is the faster one in FORTRAN
				  do k = kk, kk+N1-1
					rcohs%RC(i,k,j) = prc(k-kk+1,j-kk+1)
					if (j/=k) then
					  rcohs%RC(i,k,j) = rcohs%RC(i,k,j)*exp &
					  (-(all_goft(current_s_block%gindex(k-kk+1))%gg(gt(1)*(i-1)+1)+  &
					  conjg(all_goft(current_s_block%gindex(j-kk+1))%gg(gt(1)*(i-1)+1)))&
					  -(0.0_dp,1.0_dp)*(current_s_block%en(k-kk+1)-current_s_block%en(j-kk+1))*t)
					end if
					evops(kk,kk)%Uee(k-kk+1,j-kk+1,a,b,i) = rcohs%RC(i,k,j)
				  end do
				end do

				if (i == Nt(1)) exit

				call make_step_popcoh(t,prc,dprc)

			end do

		  end do
		end do

			deallocate(prc,dprc,cpl,RelMa)
			kk = kk + N1

			if (.not.resources_have_next_block()) exit

			call resources_next_block()

		end do

	end subroutine weak_excitons_pop_coh


	!
	! Calculates evolution of 1-to-2 exciton state coherences
	!
	subroutine weak_excitons_fe_coh()

		integer  :: kk,i,a,c,B,D,s,p,l,m,N2
		real(dp) :: t
		complex(dpc), dimension(:,:,:), allocatable :: fc, dfc

		call resources_rewind_blocks()

		kk = 1
		! loop over blocks
		do

			N2 = N1*(N1-1)/2				! to be moved into resources.f90 later
			allocate(fc(N1,N1,N1),dfc(N1,N1,N1))

			allocate(cpl(Nt(1),N1,N1)) ! preparing matrix of couplings
			cpl = 0.0_dp
			call create_couplings()
			allocate(RelMa(N1,N1,N1,N1,Nt(1))) ! preparing relaxation matrix
			RelMa = 0.0_dp
			call create_relaxation()

			evops(kk,kk)%Ufe = 0.0_dp

			! loop over i.c.'s
			do a = 1, N1
			   B = 0
			 do s = 1,   N1
			 do p = s+1, N1
				B = B+1

				fc = 0.0_dp; dfc = 0.0_dp
				fc(a,s,p) = 1.0_dp

				i = 0
				! loop over time
				do

				  i = i + 1
				  t = (i-1)*gt(1)*dt

				  ! save coherences re-indexing two-exiton part
				  do c = 1, N1
					 D = 0
				  do l = 1,   N1
				  do m = l+1, N1
					 D = D+1

!					evops(kk,kk)%Ufe(c,D,a,B,i) = fc(c,l,m) * &
					if (c == l) then
						evops(kk,kk)%Ufe(D,c,B,a,i) = fc(c,l,m) * &
						exp(  &
						  -conjg(all_goft(current_s_block%gindex(m))%gg(gt(1)*(i-1)+1))&
						  -(0.0_dp,1.0_dp)*(current_s_block%en(c)-current_s_block%en(l) &
						  -current_s_block%en(m)+rwa)*t)
					else if (c == m) then
						evops(kk,kk)%Ufe(D,c,B,a,i) = fc(c,l,m) * &
						exp(  &
						  -conjg(all_goft(current_s_block%gindex(l))%gg(gt(1)*(i-1)+1)) &
						  -(0.0_dp,1.0_dp)*(current_s_block%en(c)-current_s_block%en(l) &
						  -current_s_block%en(m)+rwa)*t)
					else
						evops(kk,kk)%Ufe(D,c,B,a,i) = fc(c,l,m) * &
						exp( -(all_goft(current_s_block%gindex(c))%gg(gt(1)*(i-1)+1)  &
						  +conjg(all_goft(current_s_block%gindex(l))%gg(gt(1)*(i-1)+1)) &
						  +conjg(all_goft(current_s_block%gindex(m))%gg(gt(1)*(i-1)+1)))&
						  -(0.0_dp,1.0_dp)*(current_s_block%en(c)-current_s_block%en(l) &
						-current_s_block%en(m)+rwa)*t)
					end if
! original implementation
!					evops(kk,kk)%Ufe(D,c,B,a,i) = fc(c,l,m) * &
!					exp( -(all_goft(current_s_block%gindex(c))%gg(gt(1)*(i-1)+1)  &
!					+conjg(all_goft(current_s_block%gindex(l))%gg(gt(1)*(i-1)+1)) &
!					+conjg(all_goft(current_s_block%gindex(m))%gg(gt(1)*(i-1)+1)))&
!					-(0.0_dp,1.0_dp)*(current_s_block%en(c)-current_s_block%en(l) &
!					-current_s_block%en(m)+rwa)*t )
				  end do
				  end do
				  end do

				  if (i == Nt(1)) exit

				  call make_step_fecoh(t,fc,dfc)

				end do

			  end do
			  end do
			end do

			deallocate(fc,dfc,cpl,RelMa)
			kk = kk + N1

			if (.not.resources_have_next_block()) exit

			call resources_next_block()

		end do

	end subroutine weak_excitons_fe_coh


!------------------------------------- make_step functions ----------------------------------------


	!
	! Step in coherence integration
	!
    subroutine make_step_coh(t,gc,dgc)
        real(dp), intent(in) :: t
        complex(dpc), intent(inout), dimension(:) :: gc
        complex(dpc), intent(inout), dimension(:) :: dgc
        integer :: i
        real(dp) :: tt

        tt = t
        do i = 1, gt(1)
            call derivs_coh(tt,gc,dgc)
            call ode_rk4(gc,dgc,tt,dt,gc,derivs_coh)
            tt = tt + dt
        end do

    end subroutine make_step_coh


	!
	! Step in population integration
	!
	subroutine make_step_pop(t,rho,drho)
		real(dp), intent(in) :: t
		real(dp), intent(inout), dimension(:) :: rho
		real(dp), intent(inout), dimension(:) :: drho
		integer :: i
		real(dp) :: tt

		tt = t
		do i = 1, gt(1)
			call derivs_pop(tt,rho,drho)
			call ode_rk4(rho,drho,tt,dt,rho,derivs_pop)
			tt = tt + dt
		end do

	end subroutine make_step_pop

	!
	! Step in population/non-optical coherences integration
	!
	subroutine make_step_popcoh(t,prc,dprc)
		real(dp), intent(in) :: t
		complex(dpc), intent(inout), dimension(:,:) :: prc
		complex(dpc), intent(inout), dimension(:,:) :: dprc
		integer :: i
		real(dp) :: tt

		tt = t
		do i = 1, gt(1)
			call derivs_popcoh(tt,prc,dprc)
			call ode_rk4(prc,dprc,tt,dt,prc,derivs_popcoh)
			tt = tt + dt
		end do
	end subroutine make_step_popcoh


	!
	! Step in Ufe coherences
	!
    subroutine make_step_fecoh(t,fe,dfe)
        real(dp), intent(in) :: t
        complex(dpc), intent(inout), dimension(:,:,:) :: fe
        complex(dpc), intent(inout), dimension(:,:,:) :: dfe
        integer :: i
        real(dp) :: tt

        tt = t
        do i = 1, gt(1)
            call derivs_fecoh(tt,fe,dfe)
            call ode_rk4(fe,dfe,tt,dt,fe,derivs_fecoh)
            tt = tt + dt
        end do

    end subroutine make_step_fecoh


!------------------------------------- RHS's of the EOM of the RDM --------------------------------


    !
	! Right hand side of the equations of motion for the populations
	!
	subroutine derivs_pop(x,y,dxdy)
		real(dp), intent(in)				 :: x ! time
		real(dpc), intent(in), dimension(:)  :: y
		real(dpc), intent(out), dimension(:) :: dxdy

		dxdy = matmul(K,y)

	end subroutine derivs_pop

	!
	! Right hand side of the equation of motion for the populations/NOC's
	!
	subroutine derivs_popcoh(x,y,dxdy)
		real(dp), intent(in) :: x
		complex(dpc), intent(in), dimension(:,:) :: y
		complex(dpc), intent(out), dimension(:,:) :: dxdy

		integer :: a, b, c, d
!not used it current scheme  complex(dpc), dimension(N1,N1) :: DM   ! matrix of dynamics without coupling

		! these are used for interpolation of g{dot}(t) aka hoft's:
		integer :: i0, i, j
		real(dp) :: x0,x1,dx
!		complex(dpc), dimension(2) :: fx, dfx ! with the view to having only two different g's
		complex(dpc), dimension(N1,N1)	   :: cplInt, dcplInt  ! coupling interpolation
		complex(dpc), dimension(N1,N1,N1,N1) :: relInt, drelInt  ! relaxation interpolation


		i0 = floor(x/(dt*gt(1))) + 1
		x0 = (i0-1)*dt*gt(1)
		x1 = i0*dt*gt(1)
		dx  = (x-x0)

!		! interpolating hoft's
!		do j = 1,nr_gofts
!		  if ((i0 + 1) >= Nt(1)) then
!			dfx(j) = 0.0_dp
!		  else
!			dfx(j) = (all_hoft(j)%gg(i0+1) - all_hoft(j)%gg(i0))/(dt*gt(1))
!		  end if
!		  fx(j)  = all_hoft(j)%gg(i0) + dfx(j)*dx
!		end do

		!interpolating coupling and relaxation
		if ((i0 + 1) >= Nt(1)) then
			dcplInt = 0.0_dp
			drelInt = 0.0_dp
		  else
			dcplInt = (cpl(i0+1,1:N1,1:N1) - cpl(i0,1:N1,1:N1))/(dt*gt(1))
			drelInt = (RelMa(:,:,:,:,i0+1) - RelMa(:,:,:,:,i0))/(dt*gt(1))
		end if
		cplInt = cpl(i0,1:N1,1:N1) + dcplInt*dx
		relInt = RelMa(:,:,:,:,i0) + drelInt*dx

!		do j  = 1, N1
!		  do i = 1, N1
!			  a = current_s_block%gindex(i)
!			  b = current_s_block%gindex(j)
!			if (i == j) then
!			  DM(i,j) = 0.0_dp  !no population dynamics at this step
!			else
!			  DM(i,j) = -(0.0_dp,1.0_dp)*(en(i)-en(j)) - (fx(a) + conjg(fx(b)))
!			end if
!		  end do
!		end do

		dxdy = 0.0_dp

		dxdy = -(0.0_dp,1.0_dp)*(matmul(cplInt,y) - matmul(y,cplInt));


		do d = 1, N1
		  do c = 1, N1
			do b = 1, N1
			  do a = 1, N1

				dxdy(a,b) = dxdy(a,b) - relInt(a,c,c,d)*y(d,b)		&
									  + conjg(relInt(c,a,b,d))*y(c,d) &
									  + relInt(d,b,a,c)*y(c,d)		&
									  - conjg(relInt(b,d,d,c))*y(a,c)

			  end do
			end do
		  end do
		end do

! dimeric input by hand for testing:

!		dxdy(1,1) = dxdy(1,1) - (relInt(1,2,2,1)+conjg(relInt(1,2,2,1)))*y(1,1) + &
!								(relInt(2,1,2,1)+conjg(relInt(2,1,2,1)))*y(2,2)
!		dxdy(1,2) = dxdy(1,2) - (relInt(1,2,2,1)+conjg(relInt(2,1,1,2)))*y(1,2) + &
!								(relInt(2,1,1,2)+conjg(relInt(1,2,2,1)))*y(2,1)
!		dxdy(2,1) = dxdy(2,1) - (relInt(2,1,1,2)+conjg(relInt(1,2,2,1)))*y(2,1) + &
!								(relInt(1,2,2,1)+conjg(relInt(2,1,1,2)))*y(1,2)
!		dxdy(2,2) = dxdy(2,2) - (relInt(2,1,1,2)+conjg(relInt(2,1,1,2)))*y(2,2) + &
!								(relInt(1,2,1,2)+conjg(relInt(1,2,1,2)))*y(1,1)

! testing --- --- ---
!	 open(unit=12,file=trim(file_join(out_dir,"nonoptrelax.dat")),position='append')
!	 write(12, '(3f16.8, /)') x, real(relInt(a,b,c,d)), aimag(relInt(a,b,c,d))
!				  aimag(-relInt(1,2,2,1)*(y(1,1)-y(2,1)-y(2,2))-conjg(relInt(1,2,2,1))*y(1,1)+ &
!				  conjg(relInt(2,1,2,1))*y(2,2)), "---", &
!				  real(-relInt(1,2,2,1)*(y(1,1)-y(2,1)-y(2,2))-conjg(relInt(1,2,2,1))*y(1,1)+ conjg(relInt(2,1,2,1))*y(2,2))
!	 close(unit=12)
! testing --- --- ---

	end subroutine derivs_popcoh


	!
	! R.h.s. of the EOM for fe coherences
	!
    subroutine derivs_fecoh(x,y,dxdy)
        real(dp), intent(in) :: x
        complex(dpc), intent(in), dimension(:,:,:) :: y
        complex(dpc), intent(out), dimension(:,:,:) :: dxdy
        integer :: a, c, s, p
        integer :: i0, i, j
        real(dp) :: x0,x1,dx
        complex(dpc), dimension(N1,N1)       :: cplInt, dcplInt  ! coupling interpolation
        complex(dpc), dimension(N1,N1,N1,N1) :: relInt, drelInt  ! relaxation interpolation
        complex(dpc), dimension(N1) :: pom

        i0 = floor(x/(dt*gt(1))) + 1
        x0 = (i0-1)*dt*gt(1)
        x1 = i0*dt*gt(1)
        dx = (x-x0)

        !interpolating coupling and relaxation
        if ((i0 + 1) >= Nt(1)) then
            dcplInt = 0.0_dp
            drelInt = 0.0_dp
          else
            dcplInt = (cpl(i0+1,1:N1,1:N1) - cpl(i0,1:N1,1:N1))/(dt*gt(1))
            drelInt = (RelMa(:,:,:,:,i0+1) - RelMa(:,:,:,:,i0))/(dt*gt(1))
        end if
        cplInt = cpl(i0,1:N1,1:N1) + dcplInt*dx
        relInt = RelMa(:,:,:,:,i0) + drelInt*dx

        dxdy = 0.0_dp

        do s = 1,   N1
        do p = s+1, N1
        do a = 1,   N1

          dxdy(a,s,p) = (0.0_dp,1.0_dp)*(         &
            - sum(cplInt(a,:)*y(:,s,p))            &
            + sum(y(a,s,s+1:N1)*cplInt(s+1:N1,p)) &
            + sum(y(a,1:p-1,p)*cplInt(1:p-1,s))   &
            + sum(y(a,p,p+1:N1)*cplInt(p+1:N1,s)) &
            + sum(y(a,1:s-1,s)*cplInt(1:s-1,p))  )!)

          do c = 1, N1
                                ! SECOND ORDER PART
          dxdy(a,s,p) = dxdy(a,s,p) - sum(relInt(a,c,c,:)*y(:,s,p))  + (                     &
            sum((conjg(relInt(c,a,p,1:s-1 )) + relInt(1:s-1 ,p,a,c))*y(c,1:s-1,s     ))    + &
            sum((conjg(relInt(c,a,p,s+1:N1)) + relInt(s+1:N1,p,a,c))*y(c,s    ,s+1:N1))    + &
            sum((conjg(relInt(c,a,s,p+1:N1)) + relInt(p+1:N1,s,a,c))*y(c,p,p+1:N1    ))    + &
            sum((conjg(relInt(c,a,s,1:p-1 )) + relInt(1:p-1 ,s,a,c))*y(c,1:p-1,p     ))      )

            ! c = alpha
          dxdy(a,s,p) = dxdy(a,s,p) -  (                     &
            2.0_dp*sum(conjg(relInt(p,c+1:N1,s,c)+relInt(s,c,p,c+1:N1)+relInt(p,c,s,c+1:N1)+relInt(s,c+1:N1,p,c))*y(a,c,c+1:N1))      &
           )

          end do

          ! (1)
          pom = 0.0_dp
          do i = s+1,N1
          pom(s+1:N1) = pom(s+1:N1) + conjg(relInt(p,i,i,s+1:N1))
          end do
          dxdy(a,s,p) = dxdy(a,s,p) - sum(pom(s+1:N1)*y(a,s,s+1:N1))

          ! (2)
          pom = 0.0_dp
          do i = s+1,N1
          pom(1:s-1) = pom(1:s-1) + conjg(relInt(p,i,i,1:s-1))
          end do
          dxdy(a,s,p) = dxdy(a,s,p) - sum(pom(1:s-1)*y(a,1:s-1,s))

          ! (3)
          pom = 0.0_dp
          do i = 1,p-1
          pom(p+1:N1) = pom(p+1:N1) + conjg(relInt(s,i,i,p+1:N1))
          end do
          dxdy(a,s,p) = dxdy(a,s,p) - sum(pom(p+1:N1)*y(a,p,p+1:N1))

          ! (4)
          pom = 0.0_dp
          do i = 1,p-1
          pom(1:p-1) = pom(1:p-1) + conjg(relInt(s,i,i,1:p-1))
          end do
          dxdy(a,s,p) = dxdy(a,s,p) - sum(pom(1:p-1)*y(a,1:p-1,p))

          ! (5)
          pom = 0.0_dp
          do i = 1,s-1
          pom(1:s-1) = pom(1:s-1) + conjg(relInt(p,i,i,1:s-1))
          end do
          dxdy(a,s,p) = dxdy(a,s,p) - sum(pom(1:s-1)*y(a,1:s-1,s))

          ! (6)
          pom = 0.0_dp
          do i = 1,s-1
          pom(s+1:N1) = pom(s+1:N1) + conjg(relInt(p,i,i,s+1:N1))
          end do
          dxdy(a,s,p) = dxdy(a,s,p) - sum(pom(s+1:N1)*y(a,s,s+1:N1))

          ! (7)
          pom = 0.0_dp
          do i = p+1,N1
          pom(p+1:N1) = pom(p+1:N1) + conjg(relInt(s,i,i,p+1:N1))
          end do
          dxdy(a,s,p) = dxdy(a,s,p) - sum(pom(p+1:N1)*y(a,p,p+1:N1))

          ! (8)
          pom = 0.0_dp
          do i = p+1,N1
          pom(1:p-1) = pom(1:p-1) + conjg(relInt(s,i,i,1:p-1))
          end do
          dxdy(a,s,p) = dxdy(a,s,p) - sum(pom(1:p-1)*y(a,1:p-1,p))


        end do
        end do
        end do

    end subroutine derivs_fecoh


	!
	! Right hand side of the equation of motion for the coherences
	!
	subroutine derivs_coh(x,y,dxdy)
		real(dp), intent(in)					:: x ! time
		complex(dpc), intent(in), dimension(:)  :: y
		complex(dpc), intent(out), dimension(:) :: dxdy

		integer :: a,c	! what are these for?.. (seem to be redundant)
		integer :: i0, i, j
		real(dp) :: x0,x1,dx
		complex(dpc), dimension(2) :: fx, dfx
		complex(dpc), dimension(N1,N1)	   :: cplInt, dcplInt  ! coupling interpolation
		complex(dpc), dimension(N1,N1,N1,N1) :: relInt, drelInt  ! relaxation interpolation


		i0 = floor(x/(dt*gt(1))) + 1
		x0 = (i0-1)*dt*gt(1)
		x1 = i0*dt*gt(1)
		dx = (x-x0)

		!print *, "nr_gofts: ", nr_gofts

!!! we don't use these in the current scheme - we use the Relaxation Matrix instead (look further)
!		do j = 1,nr_gofts
!		  !print *, nr_gofts, j, x
!		  if ((i0 + 1) >= Nt(1)) then
!			dfx(j) = 0.0_dp
!		  else
!			dfx(j) = (all_hoft(j)%gg(i0+1) - all_hoft(j)%gg(i0))/(dt*gt(1))
!		  end if
!		  fx(j)  = all_hoft(j)%gg(i0) + dfx(j)*dx
!		end do

		!interpolating coupling and relaxation
		if ((i0 + 1) >= Nt(1)) then
			dcplInt = 0.0_dp
			drelInt = 0.0_dp
		  else
			dcplInt = (cpl(i0+1,1:N1,1:N1) - cpl(i0,1:N1,1:N1))/(dt*gt(1))
			drelInt = (RelMa(:,:,:,:,i0+1) - RelMa(:,:,:,:,i0))/(dt*gt(1))
		end if
		cplInt = cpl(i0,1:N1,1:N1) + dcplInt*dx
		relInt = RelMa(:,:,:,:,i0) + drelInt*dx

!!! IMPORTANT: by aplying time-dependent coupling, the RHS is changed as well !!!
		dxdy = 0.0_dp
		dxdy =								  &	 !-(0.0_dp,1.0_dp)*(en-rwa)*y &
			   -(0.0_dp,1.0_dp)*matmul(cplInt,y)	  ! we changed jj -> cplInt in matmul
!
! Temporarily commented out
!
        do i = 1, N1
          do j = 1, N1
     !       if (j /= i) then
                                ! SECOND ORDER PART
              dxdy(i) = dxdy(i) - sum(relInt(i,j,j,:)*y) !- relInt(i,j,j,j)*y(j))

     !         print *, relInt(i,j,j,j)
     !       end if
          end do
        end do
!stop
!		dxdy(1) = dxdy(1) - relInt(1,2,2,1)*y(1)
!		dxdy(2) = dxdy(2) - relInt(2,1,1,2)*y(2)

!		! adding the damping from bath part of Liouvillian within the previous scheme
!		do i = 1, N1
!			j = current_s_block%gindex(i)
!			dxdy(i) = dxdy(i) - fx(j)*y(i) ! 0.003_dp*y(i)  ! here g_dot_of_t
!		end do

	end subroutine derivs_coh


!------------------------------------- "auxilary" functions ---------------------------------------



	!************************************************************************
	! Creates the matrix of rates for the master equation for the populations
	!************************************************************************
	subroutine create_rates()
		real(dp) :: fwhm, W, C, test ! variables for the Gaussian correlation function
		real(dp) :: gamma
		integer :: i, j

		allocate(K(N1,N1))

		gamma = 2.0_dp/300	 ! 2 over 300fs so that at w=150cm^(-1) we have downward rate of 1/300
		fwhm  = Energy_cm_to_internal * 150

		do i = 1, N1		   ! enumerate rows
		  do j = i + 1, N1	 !		   columns
			  W = current_e_block%en(i) - current_e_block%en(j)
			  call gaussn(0.0_dp,fwhm,W,C)
			  K(i,j) = gamma * C
			  if (W>0) then
				  K(i,j) = K(i,j)*exp(-W/(temp*kB_intK)) ! W>0 here!
			  end if
			  K(j,i) = K(i,j)*exp(W/(temp*kB_intK)) ! if W>0 - compensates, if W<0 - detailed ballance
		  end do
		end do

		do i = 1, N1		   ! now enumerate columns
		  K(i,i) = 0		   !			   rows
		  do j = 1, N1
			  if (i .ne. j) then
				K(i,i) = K(i,i) - K(j,i)
			  end if
		  end do
		end do
!testing part  ---  requires removal at later time:
!		open(unit=12,file=trim(file_join(out_dir,"rates.dat")))

!		W = current_e_block%en(1) - current_e_block%en(2)
!		call gaussn(0.0_dp,fwhm,W,C)
!		test = exp(W/(temp*kB_intK)) ! here W<0; we know from data
!		write(12,'(f12.8, 1x, f12.8, 1x, f12.8)') K(1,1), K(1,2), K(1,3)
!		write(12,'(f12.8, 1x, f12.8, 1x, f12.8)') K(2,1), K(2,2), K(2,3)
!		write(12,'(f12.8, 1x, f12.8, 1x, f12.8)') K(3,1), K(3,2), K(3,3)

!		close(unit=12)
!testing part  ---
	end subroutine create_rates



	!**********************************************************************
	! Creates a matrix of time-dependent couplings between sites in a block
	!**********************************************************************
	subroutine create_couplings()
		integer :: i, j, k  ! intention: loops
		integer :: a, b	    ! intention: goft indices
		real(dp) :: t

		do k  = 1, N1
		  do j = 1, N1
			  a = current_s_block%gindex(j)
			  b = current_s_block%gindex(k)
			  if (j == k) then
				cpl(1:Nt(1),j,k) = 0.0_dp
			  else
				do i = 1, Nt(1)
					t = (i-1)*gt(1)*dt  ! if 'physical' time is needed in the calculation, put into 'i' loop
					cpl(i,j,k) = jj(j,k) * exp((0.0_dp,1.0_dp)*(en(j)-en(k))*t - &
					             conjg(all_goft(a)%gg(gt(1)*(i-1)+1))-all_goft(b)%gg(gt(1)*(i-1)+1))
				end do
			  end if
		  end do
		end do
! testing --- --- ---
!		open(unit=12,file=trim(file_join(out_dir,"goft"))) ! coupling12.dat
!		do i = 1, Nt(1)
!		  write(12, '(f16.8, 1x, f16.8)') real(all_goft(1)%gg(i)), aimag(all_goft(1)%gg(i))
!		  !real(cpl(i,1,2))! i*gt(1)*dt, aimag(cpl(i,1,2))
!		end do
!		close(unit=12)
!		do i = 1, 5
!		t = (i-1)*gt(1)*dt
!		print *, "cpl(1,2)/jj(1,2) in cpl: ", &
!		            exp((0.0_dp,1.0_dp)*(en(1)-en(2))*t-conjg(all_goft(1)%gg(i))-all_goft(1)%gg(i))
!		end do
! testing --- --- ---

	end subroutine create_couplings



	!***********************************************
	! Creates the ralaxation matrix RelMa(a,b,c,d,t)
	!***********************************************
	subroutine create_relaxation()

		integer :: i

		allocate(R1(N1,N1,N1,N1,Nt(1)),R2(N1,N1,N1,N1,Nt(1)))
		R1 = 0.0_dp; R2 = 0.0_dp

		call relax_delta() !relax_any() !

		RelMa = R1 + R2

	!!! TEST !!!
	open(unit=11,file=trim(file_join(out_dir,"RelMa.dat")))
		do i = 1, Nt(1)
		  write(11,'(3f16.8)') gt(1)*(i-1)*dt, real(RelMa(1,2,1,2,i)), aimag(RelMa(1,2,1,2,i))
!		  									   real(R1(1,2,1,2,i)), aimag(R1(1,2,1,2,i))!, &
!		  									   real(R2(1,2,2,1,i)), aimag(R2(1,2,2,1,i))
		end do
	close(unit=11)
	!!! TEST !!!

		deallocate(R1,R2)

	end subroutine create_relaxation


	!
	! if in the SSF MODE_TYPE = DELTA
	!
	subroutine relax_delta()
		integer :: i, a, b, c, d  ! intention: loops
		real(dp) :: t			  ! intention: 'physical' time
		complex(dpc) :: aux
		real(dp), dimension(N1,N1,N1,N1) :: coefA, coefB, coefC !alpha,beta,gamma of the homog.lim.
		real(dp), dimension(N1,N1) :: W  ! matrix of frequencies W(a,b)
		real(dp), dimension(N1) :: gamma ! goft = gamma * t
		integer(dp), dimension(N1,N1) :: delta  ! Kronecker's delta; should be generalized somehow

		delta = 0
		forall(i = 1:N1)
			delta(i,i) = 1
		end forall

		! reading gamma's
		do i = 1, N1
		  a = current_s_block%gindex(i)
		  gamma(i) = real(all_hoft(a)%gg(1)) ! gg(t)=const for all times, so I take 1st
		end do

		do d = 1, N1
		  do c = 1, N1
			W(c,d) = current_s_block%en(c) - current_s_block%en(d)
			do b = 1, N1
			  do a = 1, N1
				coefA(a,b,c,d) = (1.0_dp + delta(a,c) - delta(a,d))*gamma(a) +		 &
								 (1.0_dp - delta(b,c) + delta(b,d))*gamma(b)
				coefB(a,b,c,d) = (delta(a,d) - delta(a,c))*gamma(a) + (delta(b,c) - delta(b,d))*gamma(b)
				coefC(a,b,c,d) = gamma(c) + gamma(d) - (delta(a,d) - delta(a,c))*gamma(a) -				 &
								 (delta(b,c) - delta(b,d))*gamma(b)
			  end do
			end do
		  end do
		end do

		do i = 1, Nt(1)
		  t = (i-1)*gt(1)*dt ! 'physical time'
		  do d = 1, N1
			do b = 1, N1
			  do c = 1, N1
				do a = 1, N1

				  if ((a/=b).and.(c/=d)) then
				  R2(a,b,c,d,i) = jj(a,b)*jj(c,d)*exp(-(gamma(a)+gamma(b)-(0.0_dp,1.0_dp)*W(a,b))*t)*	  &
								  (1.0_dp/(gamma(c)+gamma(d)-(0.0_dp,1.0_dp)*W(c,d)))*					  &
								  (exp(-(gamma(c)+gamma(d)-(0.0_dp,1.0_dp)*W(c,d))*t)-1.0_dp)

				  if (coefC(a,b,c,d)-coefB(a,b,c,d)-(0.0_dp,1.0_dp)*W(c,d) == 0.0_dp) then
						aux = t
				  else
						aux = (1.0_dp/(coefC(a,b,c,d)-coefB(a,b,c,d)-(0.0_dp,1.0_dp)*W(c,d)))*		  &
								  (exp((coefC(a,b,c,d)-coefB(a,b,c,d)-(0.0_dp,1.0_dp)*W(c,d))*t)-1.0_dp)
				  end if
				  R1(a,b,c,d,i) = jj(a,b)*jj(c,d)*														    &
								  exp((0.0_dp,1.0_dp)*(W(a,b)+W(c,d))*t-(coefA(a,b,c,d)+coefC(a,b,c,d))*t)* &
								  aux
				  end if

				end do
			  end do
			end do
		  end do
		end do

	end subroutine relax_delta


	!
	! if in the SSF MODE_TYPE = any type of correlation
	!
	subroutine relax_any()
		complex(dpc), dimension(N1,N1,N1,N1,1024) :: M1, M2, M3 !1024 <-> Nt(1)
!		complex(dpc), dimension(1,2*Nt(1))   :: conv_out ! convolution result = 1 function
!		complex(dpc), dimension(2,2*Nt(1))   :: conv_in  ! convolution input = 2 functions
		complex(dpc), dimension(Nt(1),1)   :: ReIntCpl, ImIntCpl
		complex(dpc), dimension(Nt(1),2)   :: test_in, test_out ! testing 'behaviour' of convolution
		real(dp), dimension(Nt(1))         :: time
		real(dp), dimension(1)             :: w
		real(dp)						   :: G ! in long times g(t) -> G*t - i*lambda*t
		integer :: i, a, b, c, d    ! intention: loops
		integer :: j, k, l, m	    ! intention: goft indices
		integer :: n_conv			! limit of extent for convolution; BEWARE!
		real(dp), dimension(N1,N1) :: delta  ! Kronecker's delta; should be generalized somehow !!! 4 <-> N1 after testing

		integer(dp) :: d_ratio ! ratio of internal time step dt to internal time step dt1 during <densification>
		complex(dpc), dimension(:),   allocatable :: M2e, M3e
		complex(dpc), dimension(:,:), allocatable :: conv_in_e, conv_out_e

		! n_conv=1024 <-> Nt(1)
		n_conv = 1024
		if (Nt(1) < n_conv) then
			print *, "ERROR: Insufficient time span: Nt(1) <n_conv!"
			stop
		end if

		d_ratio = 64 ! 64; must be a power of 2
		allocate(M2e(n_conv*d_ratio),M3e(n_conv*d_ratio))
		! (Nt(1);2*Nt(1)] FT 'sees' it as t < 0 => leave zeros there
		allocate(conv_in_e(2,2*n_conv*d_ratio),conv_out_e(1,2*n_conv*d_ratio))

		delta = 0.0_dp
		forall(i = 1:N1)
		  delta(i,i) = 1.0_dp
		end forall

		forall(i = 1:Nt(1))
		  time(i) = (i-1)*gt(1)*dt
		end forall

		M1 = 0.0_dp; M2 = 0.0_dp; M3 = 0.0_dp
		w = 0.0_dp

		do i = 1, n_conv !1024 <-> Nt(1)
			do d = 1, N1
			do c = 1, N1
			if (c/=d) then
				do b = 1, N1
				do a = 1, N1
				if (a/=b) then

				  j = current_s_block%gindex(a)
				  k = current_s_block%gindex(b)
				  l = current_s_block%gindex(c)
				  m = current_s_block%gindex(d)

				  M1(a,b,c,d,i) = exp(-conjg(all_goft(j)%gg(gt(1)*(i-1)+1)) - all_goft(k)%gg(gt(1)*(i-1)+1)- &
				  				  delta(a,c)*all_goft(j)%gg(gt(1)*(i-1)+1) + &
				  				  delta(a,d)*all_goft(j)%gg(gt(1)*(i-1)+1) + &
				  				  delta(b,c)*all_goft(k)%gg(gt(1)*(i-1)+1) - &
				  				  delta(b,d)*all_goft(k)%gg(gt(1)*(i-1)+1))
				  M3(a,b,c,d,i) = exp(-conjg(all_goft(l)%gg(gt(1)*(i-1)+1)) - all_goft(m)%gg(gt(1)*(i-1)+1)- &
				  				  delta(a,c)*conjg(all_goft(j)%gg(gt(1)*(i-1)+1)) + &
								  delta(a,d)*conjg(all_goft(j)%gg(gt(1)*(i-1)+1)) + &
				  				  delta(b,c)*conjg(all_goft(k)%gg(gt(1)*(i-1)+1)) - &
				  				  delta(b,d)*conjg(all_goft(k)%gg(gt(1)*(i-1)+1)))

				  G = 0.0_dp
				  if ((d==b).AND.(c==a)) then
						G = real( all_goft(j)%gg(gt(1)*Nt(1)) - all_goft(j)%gg(gt(1)*(Nt(1)-10)) + &
								  all_goft(k)%gg(gt(1)*Nt(1)) - all_goft(k)%gg(gt(1)*(Nt(1)-10)) ) &
							/( gt(1)*Nt(1) - gt(1)*(Nt(1)-10) )
						M2(a,b,c,d,i) = exp(all_goft(j)%gg(gt(1)*(i-1)+1) + all_goft(k)%gg(gt(1)*(i-1)+1) &
										 - G*time(i))
				  elseif (c==a) then
						G = real( all_goft(j)%gg(gt(1)*Nt(1)) - all_goft(j)%gg(gt(1)*(Nt(1)-10)) )&
							/( gt(1)*Nt(1) - gt(1)*(Nt(1)-10) )
						M2(a,b,c,d,i) = exp(all_goft(j)%gg(gt(1)*(i-1)+1) - G*time(i))
				  elseif (d==b) then
						G = real( all_goft(k)%gg(gt(1)*Nt(1)) - all_goft(k)%gg(gt(1)*(Nt(1)-10)) )&
							/( gt(1)*Nt(1) - gt(1)*(Nt(1)-10) )
						M2(a,b,c,d,i) = exp(all_goft(k)%gg(gt(1)*(i-1)+1) - G*time(i))
				  else
						M2(a,b,c,d,i) = exp(delta(a,c)*all_goft(j)%gg(gt(1)*(i-1)+1) - &
											delta(a,d)*all_goft(j)%gg(gt(1)*(i-1)+1) - &
											delta(b,c)*all_goft(k)%gg(gt(1)*(i-1)+1) + &
											delta(b,d)*all_goft(k)%gg(gt(1)*(i-1)+1))
				  end if

				  M3(a,b,c,d,i) = M3(a,b,c,d,i)*exp(  (-G + (0.0_dp,1.0_dp)*(en(c)-en(d)) ) * time(i)  )
				  ! Safety trigger for positive exponential: EXP[150]
				  if ( (G*time(i)) < 150.0_dp ) then
				    M1(a,b,c,d,i) = M1(a,b,c,d,i)*exp(  ( G + (0.0_dp,1.0_dp)*(en(a)-en(b)) ) * time(i)  )
				  else
					M1(a,b,c,d,i) = (0.0_dp,0.0_dp)
				  end if

				end if
				end do
				end do
			end if
			end do
			end do
		end do

		do d = 1, N1
  		do c = 1, N1
		if (c/=d) then

			call primitive(time, real(cpl(1:Nt(1),c,d)),w,ReIntCpl(1:Nt(1),1))
			call primitive(time,aimag(cpl(1:Nt(1),c,d)),w,ImIntCpl(1:Nt(1),1))

			do b = 1, N1
			do a = 1, N1
			if (a/=b) then

				  R2(a,b,c,d,1:Nt(1)) = -cpl(1:Nt(1),a,b) * &
				  						(ReIntCpl(1:Nt(1),1) + (0.0_dp,1.0_dp)*ImIntCpl(1:Nt(1),1))

!----------------------------------------------- <densification> of M2, M3 prior to convolution
				  call densify_linear(M2(a,b,c,d,1:n_conv),M2e)
				  call densify_linear(M3(a,b,c,d,1:n_conv),M3e)
				  conv_in_e = 0.0_dp; conv_out_e = 0.0_dp
!-----------------------------------------------
				  conv_in_e(1,1:n_conv*d_ratio) = M2e
				  conv_in_e(2,1:n_conv*d_ratio) = M3e
				  call fft_row_convolution(conv_out_e(1:1,1:2*n_conv*d_ratio), &
				  						   conv_in_e(1:1,1:2*n_conv*d_ratio),conv_in_e(2:2,1:2*n_conv*d_ratio), &
										   NOSE_QME_NZ_REG_TIME_PER_TMAX)
!----------------------------------------------- picking out values of covolution rezult
!				  R1(a,b,c,d,1:Nt(1)) =  jj(a,b)*jj(c,d)*M1(a,b,c,d,1:Nt(1))*conv_out(1,1:Nt(1))*gt(1)*dt ! t.step
				  do i = 1, n_conv ! Nt(1) <-> n_conv
					R1(a,b,c,d,i) =  jj(a,b)*jj(c,d)*M1(a,b,c,d,i)*conv_out_e(1,(i-1)*d_ratio+1)*(gt(1)*dt/d_ratio)
				  end do
				  R1(a,b,c,d,n_conv+1:Nt(1)) = R1(a,b,c,d,n_conv)
!-----------------------------------------------
		 	end if
			end do
			end do

		end if
		end do
		end do

!*************************************      testing      ******************************************************
		open(unit=12,file=trim(file_join(out_dir,"conv.dat")))
!		call fft_row_regularized(M2(1:1,1,1,1,:),1,NOSE_QME_NZ_REG_TIME_PER_TMAX)
!		call fft_row_convolution(conv(1:1,1:Nt(1)),M2(1:1,2,2,1,1:Nt(1)),M3(1:1,2,2,1,1:Nt(1)), &
!		                         NOSE_QME_NZ_REG_TIME_PER_TMAX)
!		G = real( all_goft(1)%gg(gt(1)*Nt(1)) - all_goft(1)%gg(gt(1)*(Nt(1)-10)) ) &
!		/( gt(1)*Nt(1) - gt(1)*(Nt(1)-10) )
!				  conv_in(1,1:Nt(1)) = M2(1,2,1,2,1:Nt(1))
!				  conv_in(2,1:Nt(1)) = M3(1,2,1,2,1:Nt(1))
!				  call fft_row_convolution(conv_out(1:1,1:2*Nt(1)),conv_in(1:1,1:2*Nt(1)),conv_in(2:2,1:2*Nt(1)), &
!										   NOSE_QME_NZ_REG_TIME_PER_TMAX)
		do i = 1, n_conv ! Nt(1)
		  write(12, '(3f18.8)') (i-1)*dt*gt(1), &
									real(R1(1,2,2,1,i)), real(R1(2,1,1,2,i))
!		 		real(M1(1,2,2,1,i)), aimag(M1(1,2,2,1,i)), &
!				real(M2(1,2,2,1,i)), aimag(M2(1,2,2,1,i)) , &
!				real(M3(1,2,2,1,i)), aimag(M3(1,2,2,1,i)) !, &
!		  write(12, '(4f18.8)') (i-1)*dt*gt(1), real(all_goft(1)%gg(gt(1)*(i-1)+1)), &
!		  						aimag(all_goft(1)%gg(gt(1)*(i-1)+1)), G*(i-1)*dt*gt(1)
!		  real(test_out(i,1)), real(test_out(i,2)) !, real(M2(1,1,1,1,i)), aimag(M2(1,1,1,1,i))
		end do
		close(unit=12)
!***********************************************************************************************************

		deallocate(M2e, M3e, conv_in_e, conv_out_e)

	end subroutine relax_any

	!
	! Densifying a function by d_ratio using simple linear interpolation
	!
	subroutine densify_linear(in,out)
 		complex(dpc), dimension(:), intent(in)		:: in
 		complex(dpc), dimension(:), intent(out)		:: out
 		integer :: i, j, ratio, size_in, size_out

 		size_in  = size(in,1)
 		size_out = size(out,1)
 		ratio = size_out/size_in

 		do i = 1, size_in
 		  do j = 1, ratio
 		  	if (i/=size_in) then
 		      out((i-1)*ratio+j) = in(i) + (in(i+1)-in(i))*(j-1)/ratio
 		      else
 		      out((i-1)*ratio+j) = in(i)
 		    end if
 		  end do
 		end do

	end subroutine


end module qme_weak_excitons
