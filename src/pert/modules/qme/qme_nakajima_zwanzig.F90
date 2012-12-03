!
!
!
!
#include "util_allocation.h"

#define SUPERINDEX_FROM_K_L(k, l, lmax) ((l) + (lmax)*((k)-1))
#define K_FROM_SUPERINDEX(superindex, lmax) (((superindex) - 1) / (lmax) + 1)
#define L_FROM_SUPERINDEX(superindex, lmax) (mod(((superindex)-1) , (lmax)) + 1)

#define PALLOCATABLE pointer
!allocatable
!pointer
#define PALLOCATED   associated
!allocated
!associated
#define PNULLIFY(x) nullify(x)
! leave empty if not using pointers

!
!
!
!
module qme_nakajima_zwanzig

	use prepare
	use twod

	use resources_qme
	use std_types
	use nakajima_zwanzig_shared

	use numer_ode
	use numer_fft
	use numer_matrix
	use sci_misc

	use util_allocation

	implicit none

 	! declarations

	type exc1_characteristics ! characteristics for each 1exciton
	    complex(dpc), dimension(:,:), PALLOCATABLE 		:: K_exc => NULL()
		complex(dpc), dimension(:,:,:), PALLOCATABLE 		:: dLambda_exc => NULL()
		complex(dpc), dimension(:,:,:), PALLOCATABLE 		:: dLambda_exc_w => NULL()

		complex(dpc), dimension(:,:,:), PALLOCATABLE 		:: Lambda_exc => NULL()
	end type

	type exc2_characteristics ! characteristics for each 2exciton
		complex(dpc), dimension(:,:), PALLOCATABLE 		:: K_2_exc => NULL()
		complex(dpc), dimension(:,:,:), PALLOCATABLE 		:: dLambda_2_exc => NULL()
		complex(dpc), dimension(:,:,:), PALLOCATABLE 		:: dLambda_2_exc_w => NULL()

		complex(dpc), dimension(:,:,:), PALLOCATABLE 		:: Lambda_2_exc => NULL()
	end type

	type(exc1_characteristics), private, dimension(:), allocatable	:: exc1
	type(exc2_characteristics), private, dimension(:), allocatable	:: exc2

	! Fourier transform of heaviside theta to ensure M(t) = 0 for t<0.
	! Frequency independent members of M have to have heaviside_w w-dependence.
	complex(dpc), private, dimension(:,:), allocatable :: heaviside_w

	integer(i4b), private	:: Nl1, Nl2

	! global variables for runge-kutta usage
	character, private :: global_type, global_submethod
	complex(dpc), private, dimension(:,:), allocatable :: global_superop_L

	! global variable for "phenomenological" ralaxation
	complex(dpc), private, dimension(:,:), allocatable :: phenom_Rates

!	! global variables tau-projectors
!	real(dp), private 				:: tau_of_projector

	! c_w for direct Lambda calculation
	complex(dpc), private, dimension(:), allocatable :: c_w

	logical, private :: Redfield_printed = .false.

	!general functions and functions of qme (methods 'q', 'Q')
	public::fill_evolution_superoperator_nakajima_zwanzig
	private::perform_superoperator
	private::minus_w_index
	private::minus_w_index_real
	private::interpolate
	private::prepare_Hamiltonian_exc
	private::prepare_Ue_exc
	private::prepare_K_exc
	private::prepare_K_excs
	private::prepare_dLambda_exc
	private::prepare_dLambdas
	private::fourier_dLambdas
	private::calculate_Delta_w
	private::generate_warnings
	private::fill_fourier_parameters
	private::secularize_superoperator
	private::fill_superoperator_M_w
	private::fill_superoperator_L
	private::init_qme
	private::deinit_qme
	private::perform_qme

	! Redfield functions (methods 'r', 's', 'R', 'S', 'p', 'P')
	private::fill_superoperator_R
	private::perform_redfield
	private::calculate_derivative_redf
	private::integrate_Lambdas
	private::phenom_relax ! 'phenomenological relaxation' of populations and coherences (methods 'p', 'P')

	! tau-redfield functions (methods 'u', 'U')
	public::goft_exciton
	public::dgoft_exciton
	public::dgoft_half_exciton
	public::dgoft_site
!	private::weigth_superoperator_by_exc_goft
	private::fill_superoperator_R_tau_addition
	private::calculate_derivative_red_tau
	public::write_redfield_tensor2

	! funkce pocitajici rovnou Lambda(w)
	!private::fill_c_w

 	contains

 	subroutine perform_superoperator(superop, target, out)
 		complex(dpc), dimension(:,:), intent(in)		:: superop, target
 		complex(dpc), dimension(:,:), intent(out)		:: out

 		complex(dpc), dimension(size(target))			:: target_as_vector
 		complex(dpc), dimension(size(out))			:: out_as_vector
 		integer			:: i,j,k,l

 		if (.not.(size(superop,2) == size(target)) .or. .not.(size(superop,1) == size(out)) ) then
 			call print_error_message(-1,  'dimension error in perform_superoperator')
 			stop
 		end if

 		do i = 1, size(target, 1)
 		do j = 1, size(target, 2)
 			target_as_vector(SUPERINDEX_FROM_K_L(i,j,size(target,2))) = target(i,j)
 		end do
 		end do

 		out_as_vector = matmul(superop, target_as_vector)

 		do i = 1, size(out, 1)
 		do j = 1, size(out, 2)
 			out(i,j) = out_as_vector(SUPERINDEX_FROM_K_L(i,j,size(out,2)))
 		end do
 		end do

 	end subroutine perform_superoperator

 	pure function minus_w_index(w_index) result(minus_w_ind)
 		integer, intent(in)		:: w_index
 		integer		:: minus_w_ind

 		if(w_index == 1 .or. w_index == grid_Nt) then
			minus_w_ind = w_index
		else
			minus_w_ind = grid_Nt*2 - w_index + 2
		end if
 	end function

 	pure function minus_w_index_real(w_index) result(minus_w_ind)
 		real(dp), intent(in)		:: w_index
 		real(dp)					:: minus_w_ind

 		if(int(w_index+0.5) == 1 .or. int(w_index+0.5) == grid_Nt) then
			minus_w_ind = w_index
		else
			minus_w_ind = grid_Nt*2 - w_index + 2
		end if
 	end function

 	function interpolate(field, real_index) result(x)
 		complex(dpc), dimension(:,:,:), intent(in)	:: field
 		real(dp), intent(in)			:: real_index
 		complex(dpc), dimension(size(field,1),size(field,2))	:: x

 		integer(i4b)	:: down_index, up_index
 		real(dp)		:: weight_down

 		if(size(field,3) < 1) then
 			call print_error_message(-1, "error in interpolate()")
 		end if

 		if(real_index < 0 .or. real_index > size(field,3)+1) then
			write(*,*) 'real_index=',real_index, down_index, up_index, weight_down
 			call print_error_message(-1, "index out of range in interpolate()")
 		end if


 		down_index	= int(real_index)
 		up_index	= down_index + 1
 		weight_down	= up_index-real_index

  		if(down_index == 0) then
 			down_index	= size(field,3)
 			up_index	= 1
 			weight_down	= 1.0_dp - real_index
 		else if(up_index == size(field,3) + 1 ) then
 			up_index	= 1
 			down_index  = size(field,3)
		else if(.not. (up_index <= size(field,3) .and. down_index >= 1) ) then
			call print_error_message(-1, "unknown error in interpolate()")
 		end if

 		x = field(:,:,down_index)*weight_down + field(:,:,up_index)*(1-weight_down)

 	end function

 	subroutine prepare_Hamiltonian_exc(Ham,type)
		complex(dpc), dimension(:,:), intent(out) 	:: Ham
		character, intent(in)					 		:: type
 		integer(i4b) 									:: i,j,k

 		if(.not.(size(Ham,1) == N1_from_type(type) .and. size(Ham,2) == N2_from_type(type))) then
 			call print_error_message(-1,  'Dimension error in prepare_Hamiltonian()')
 			stop
 		end if

		if(type == 'O' .or. type == '2') then
			Ham = 0.0_dp
		else if(type == 'E') then
 			do i=1, Nl1
			do j=1, Nl1
				if(i == j) then
					Ham(i,i) = iblocks(1,1)%eblock%en(i) - rwa
				else
					Ham(i,j) = 0
				endif
			end do
			end do
		else if(type == 'F') then
 			do i=1, Nl2
			do j=1, Nl2
				if(i == j) then
					Ham(i,i) = iblocks(1,1)%eblock%en_2(i) - 2*rwa
				else
					Ham(i,j) = 0
				endif
			end do
			end do
		end if

 	end subroutine prepare_Hamiltonian_exc

 	subroutine prepare_Ue_exc(to, t, type)
 		complex(dpc), dimension(:,:), intent(out) :: to
 		double precision, intent(in) 	:: t
 		character, intent(in)			:: type
 		integer	(i4b)					:: i,j

 		if(.not.(size(to,1) == N1_from_type(type) .and. size(to,2) == N2_from_type(type))) then
 			call print_error_message(-1,  'Dimension error in prepare_Hamiltonian()')
 			stop
 		end if

 		if(type == 'E') then
 			do i=1, Nl1
			do j=1, Nl1
				if(i == j) then
					to(i,i) = exp(-(0.0_dp,1.0_dp)*(iblocks(1,1)%eblock%en(i)-rwa)*t)
				else
					to(i,j) = 0
				endif
			end do
			end do

 		else if(type == 'F') then
 			do i=1, Nl2
			do j=1, Nl2
				if(i == j) then
					to(i,i) = exp(-(0.0_dp,1.0_dp)*(iblocks(1,1)%eblock%en_2(i)-2*rwa)*t)
				else
					to(i,j) = 0
				endif
			end do
			end do
 		end if

 	end subroutine prepare_Ue_exc

 	subroutine prepare_K_exc(to, ind, type)
 		complex(dpc), dimension(:,:), intent(out) :: to
 		complex(dpc), dimension(size(to,1),size(to,2)) :: K
 		character, intent(in) :: type
 		integer, intent(in) 	:: ind
 		integer					:: i,j,N_use

 		if(type == 'E') then
 			N_use = Nl1
 		else if(type == 'F') then
 			N_use = Nl2
 		else
 			call print_error_message(-1,  'wrong type in prepare_K_exc()')
 			stop
 		end if

 		if(.not.(size(to,1) == N_use .and. size(to,2) == N_use)) then
 			call print_error_message(-1,  'dimension error in prepare_K_exc()')
 			stop
 		end if

 		do i=1, N_use
		do j=1, N_use
			if(i == j .and. i == ind) then
				K(i,i) = 1
			else

				K(i,j) = 0
			endif
		end do
		end do

		to = K
		call operator_to_exc(to, type)

 	end subroutine prepare_K_exc

 	subroutine prepare_K_excs()
 		integer(i4b) :: i

		do i = 1, Nl1
			call prepare_K_exc(exc1(i)%K_exc,i,'E')
		end do
		if(use_twoexcitons) then
 			do i = 1, Nl2
 				call prepare_K_exc(exc2(i)%K_2_exc,i,'F') ! 'F' stands for 2-excitonic block
 			end do
 		end if

 	end subroutine prepare_K_excs

 	subroutine prepare_dLambda_exc(to, ind, type)
        complex(dpc), dimension(:,:,:), intent(out)     		:: to
 		complex(dpc), dimension(size(to,1),size(to,2))		:: Ue_exc, Ue1_exc
 		character, intent(in)									:: type
 		integer, intent(in) 									:: ind

 		double precision :: t
 		integer(i4b) :: i,j,k, N_use, ind1, ind2

 		if(type == 'E') then
 			N_use = Nl1
 		else if(type == 'F') then
 			N_use = Nl2
 			ind1 = iblocks(1,1)%eblock%ione1(ind)
 			ind2 = iblocks(1,1)%eblock%ione2(ind)
 		end if

 		if(.not.(size(to,1)==N_use .and. size(to,2)==N_use)) then
			call print_error_message(-1,  'dimension error in prepare_dLambda_exc')
			stop
		end if

		to = 0.0_dp

		! Lambdas are considered to be zero for t < 0 => maximal index grid_Nt
 		do i=1,grid_Nt
 			t = (i-1)*dt

 			call prepare_Ue_exc(Ue_exc, t, type)
 			Ue1_exc = conjg(transpose(Ue_exc))


			! here, monomer-monomer correlation functions can be implemented

			if(type == 'E') then
				to(:,:,i) = matmul(matmul(Ue_exc,exc1(ind)%K_exc),Ue1_exc)

				if(i == 1) then ! correlation fction has 0 imaginary part at t = 0
 					to(:,:,i) = to(:,:,i)*real(igofts(iblocks(1,1)%sblock%gindex(ind))%goft%ct(i))
 				else
 					to(:,:,i) = to(:,:,i)*igofts(iblocks(1,1)%sblock%gindex(ind))%goft%ct(i)
 				end if

			else if(type == 'F') then
				to(:,:,i) = matmul(matmul(Ue_exc,exc2(ind)%K_2_exc),Ue1_exc)

				if(i == 1) then ! correlation fction has 0 imaginary part at t = 0
					to(:,:,i) = to(:,:,i)*real(igofts(iblocks(1,1)%sblock%gindex(ind1))%goft%ct(i) + igofts(iblocks(1,1)%sblock%gindex(ind2))%goft%ct(i))
				else
					to(:,:,i) = to(:,:,i)*(igofts(iblocks(1,1)%sblock%gindex(ind1))%goft%ct(i) + igofts(iblocks(1,1)%sblock%gindex(ind2))%goft%ct(i))
				end if

			end if

 		end do

 	end subroutine prepare_dLambda_exc

 	subroutine prepare_dLambdas()
 		integer(i4b) :: i

		! 'E' block of 1-exciton states
		do i=1, Nl1
			call prepare_dLambda_exc(exc1(i)%dLambda_exc, i, 'E')
		end do

		if(use_twoexcitons) then
			! 'F' block of 2-exciton states
			do i=1, Nl2
				call prepare_dLambda_exc(exc2(i)%dLambda_2_exc, i, 'F')
			end do
		end if
 	end subroutine prepare_dLambdas

 	subroutine fourier_dLambdas()
 		integer(i4b) 		:: i

		! 'E' block of 1-exciton states
 		do i=1, Nl1
			exc1(i)%dLambda_exc_w = exc1(i)%dLambda_exc
			call fft_row_regularized(exc1(i)%dLambda_exc_w, 1, NOSE_QME_NZ_REG_TIME_PER_TMAX)
		end do

		if(use_twoexcitons) then
			! 'F' block of 2-exciton states
			do i=1, Nl2
				exc2(i)%dLambda_2_exc_w = exc2(i)%dLambda_2_exc
				call fft_row_regularized(exc2(i)%dLambda_2_exc_w, 1, NOSE_QME_NZ_REG_TIME_PER_TMAX)
			end do
		end if
 	end subroutine fourier_dLambdas

 	function calculate_Delta_w(c,d,type) result(Delta_w_index)
 		integer(i4b), intent(in)	:: c,d
 		character, intent(in)		:: type
 		real(dp) 					:: Delta_w_index

 		integer(i4b)	:: i,j
 		real(dp)		:: Delta_omega
 		real(dp)		:: Trange,wrange,dw,gg

 		if(	c < 1 .or. c > N1_from_type(type) .or. d < 1 .or. d > N2_from_type(type) 	&
			.or. (.not.(type == 'O' .or. type == 'E' .or. type == '2'))	) then
 			call print_error_message(-1,  'Dimension error in calculate_Delta_w()')
 			stop
 		end if

		call fill_fourier_parameters(Trange,wrange,dw,gg)

 		if(type == 'O') then
			Delta_omega = - iblocks(1,1)%eblock%en(c)  + rwa
 		elseif(type == 'E') then
			Delta_omega = iblocks(1,1)%eblock%en(d) - iblocks(1,1)%eblock%en(c)
 		elseif(type == '2') then
			Delta_omega = iblocks(1,1)%eblock%en(d) - iblocks(1,1)%eblock%en_2(c)  + rwa

 		end if

 		!Delta_w_index = int(Delta_omega / dw)
 		Delta_w_index = Delta_omega / dw

 	end function

 	subroutine generate_warnings(submethod)
 		character, intent(in) :: submethod

 		integer 			:: i
 		real(dp) 			:: maximum, tmp
 		character(len=128)	:: cbuff

		if(submethod == 'q' .or. submethod == 'Q') then

		maximum = 0
		do i=1, Nl1
			tmp = maxval(abs(exc1(i)%dLambda_exc_w(:,:,grid_Nt)))/maxval(abs(exc1(i)%dLambda_exc_w))
			if(tmp > maximum) then
				maximum = tmp
			end if
		end do

		if(maximum > 0.001) then

			write(cbuff,'(f6.2)'), maximum*100/2
			cbuff = 'estimated equillibrium inacuracy: '// trim(cbuff) // '% of trace - try smaller step'
			call print_warning_message(trim(cbuff),3)

		end if

 		end if

 	end subroutine generate_warnings

 	subroutine secularize_superoperator(S, type)
 		complex(dpc), dimension(:,:), intent(inout) 	:: S
 		character, intent(in)							:: type

 		integer :: i,j,  k,l

 		if((.not. size(S,1)==size(S,2)) .or. &
 			(.not. size(S,1) == N1_from_type(type)*N2_from_type(type)) ) then
 			call print_error_message(-1,'incorrect S-size in secularize_superoperator()')
 			stop
 		end if

 		do i=1,N1_from_type(type)
 		do j=1,N2_from_type(type)

 		do k=1,N1_from_type(type)
 		do l=1,N2_from_type(type)
 			if(type=='E'.and. .not.(i==j .and. k==l .or. i==k .and. j==l) .or.	&
 				(type=='O' .or. type=='2') .and. .not. (i==k .and. j==l) ) then
 				S(	SUPERINDEX_FROM_K_L(i, j, N2_from_type(type)),  &
					SUPERINDEX_FROM_K_L(k, l, N2_from_type(type)) ) &
							= 0.0_dp
 			end if
 		end do
 		end do

 		end do
 		end do

 	end subroutine secularize_superoperator

!------------------------------------- M & L --------------------------------

 	subroutine fill_superoperator_M_w(M_w,intended_w_index,type)
 		character, intent(in) :: type
        integer, intent(in) :: intended_w_index
        complex(dpc), dimension(:,:), intent(out) :: M_w

 		complex(dpc), dimension(N1_from_type(type),N2_from_type(type)) 	&
 											:: rho, reslt
 		complex(dpc), dimension(Nl1,Nl1)	:: Hm, Hm_conjg
 		complex(dpc), dimension(Nl2,Nl2)	:: Hm_2, Hm_2_conjg
 		integer(i4b) 	:: i,j,k,l,s, first_N, second_N
 		real(dp)		:: Delta_w_index_shift, w_index, min_w_index
 		logical			:: vratit_nulu = .false.

 	if(	(.not.(size(M_w,1) == N1_from_type(type)*N2_from_type(type))) .and. &
 		(.not.(size(M_w,2) == N1_from_type(type)*N2_from_type(type))) ) then
 			call print_error_message(-1,  'Dimension error in fill_superoperator_M_w()')
 			stop
 	end if

 	if(	(intended_w_index > size(exc1(1)%dLambda_exc_w,3)) .or. &
 		(intended_w_index < 1) ) then

 		call print_error_message(-1,  'Time/frequency dimension error in fill_superoperator_M_w()')
 		stop
 	end if



	first_N 	= N1_from_type(type)
	second_N	= N2_from_type(type)

	do k=1, first_N
	do l=1, second_N
		rho = 0.0_dp
		rho(k,l) = 1.0_dp

		! Superoperator M(t) is given by following equations. But because we
		! calculate convolution M(t-tau)U(t-tau)rho(tau), we have to include
		! U(t-tau) superop. in some clever way. It is diagonal
		! in excitonic picture and therefore we only have to shift frequency
		! by Delta_w = w_k - w_l for given rho_kl

		Delta_w_index_shift = calculate_Delta_w(k,l,type)

		if(	Delta_w_index_shift + intended_w_index > grid_Nt .and. 	&
			intended_w_index < grid_Nt) then

			w_index = grid_Nt ! all functions assumed to be zero there
			!vratit_nulu = .true.
		elseif(	Delta_w_index_shift + intended_w_index < grid_Nt .and.	&
			intended_w_index > grid_Nt	) then

			w_index = grid_Nt ! all functions assumed to be zero there
			!vratit_nulu = .true.
		else
			w_index = Delta_w_index_shift + intended_w_index
		end if

		if(w_index > grid_Nt*2) then
			w_index = w_index - grid_Nt*2
		endif
		if(w_index < 1) then
			w_index = w_index + grid_Nt*2
		end if

		! this should not happen
		if(w_index > grid_Nt*2 + 1 .or. w_index < 0) then
			call print_error_message(-1,  'error in fill_superoperator_M_w()')
			stop
		end if

		min_w_index = minus_w_index_real(w_index)

		if(w_index > grid_Nt*2+1 .or. min_w_index > grid_Nt*2+1) then
			write(*,*) 'w_index=',intended_w_index, w_index, min_w_index, Delta_w_index_shift
		end if


		Hm 			= 0.0_dp
		Hm_2 		= 0.0_dp

		Hm_conjg	= 0.0_dp
		Hm_2_conjg	= 0.0_dp

		! K_exc from 'F' and K_2_exc from 'E' are zero, only terms from
		! particular block contribute

		do i = 1, Nl1

			Hm 		= Hm 		&
		+ (0,-1)*matmul(exc1(i)%K_exc,interpolate(exc1(i)%dLambda_exc_w,w_index) )

			Hm_conjg= Hm_conjg 	&
		+ transpose(conjg( 		&
		  (0,-1)*matmul(exc1(i)%K_exc,interpolate(exc1(i)%dLambda_exc_w,min_w_index)) &
		  ))

		end do


		if(type == '2') then
		do i = 1, Nl2

			Hm_2 		= Hm_2			&
		+ (0,-1)*matmul(exc2(i)%K_2_exc,interpolate(exc2(i)%dLambda_2_exc_w,w_index))
			Hm_2_conjg 	= Hm_2_conjg	&
		+ transpose(conjg( 				&
		  (0,-1)*matmul(exc2(i)%K_2_exc,interpolate(exc2(i)%dLambda_2_exc_w,min_w_index)) &
		  ))

		end do
		end if




		!debug
		!Hm = 0.0_dp
		!Hm_conjg = 0.0_dp
		!Hm(1,1)		  = heaviside_w(1,w_index)
		!Hm_conjg(1,1) = conjg(heaviside_w(1,min_w_index))





		if(type == 'O') then
			reslt = -(0,1)*matmul(Hm,rho)
		else if(type == 'E') then
			reslt = -(0,1)*(matmul(Hm,rho) - matmul(rho,Hm_conjg))
!			reslt = 0.0_dp !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			! OPEN - shortcut only - this comment can be deleted

			do s=1,Nl1
				reslt = reslt 															&
					+ matmul(matmul(exc1(s)%K_exc,rho),									&
					transpose(conjg(interpolate(exc1(s)%dLambda_exc_w,min_w_index))) ) 	&
					+ matmul(matmul(interpolate(exc1(s)%dLambda_exc_w,w_index),rho),	&
					transpose(conjg(exc1(s)%K_exc)) )

!				reslt = reslt + 												&
!					matmul(rho,													&
!					transpose(conjg(exc1(s)%dLambda_exc_w(:,:,min_w_index)))) + 	&
!					matmul(exc1(s)%dLambda_exc_w(:,:,w_index),rho)
!
!					transpose(conjg(exc1(s)%dLambda_exc_w(:,:,min_w_index)))
!					exc1(s)%dLambda_exc_w(:,:,w_index)
!
!
!					matmul(matmul(exc1(s)%K_exc,rho),							&
!					transpose(conjg(exc1(s)%dLambda_exc_w(:,:,min_w_index))))
!					matmul(matmul(exc1(s)%dLambda_exc_w(:,:,w_index),rho),		&
!					transpose(conjg(exc1(s)%K_exc)))


			end do

!				if(k == 2 .and. l == 1) then
!					reslt(1,1) = heaviside_w(1,w_index)
!					reslt(1,2) = heaviside_w(1,w_index)
!					reslt(2,1) = heaviside_w(1,w_index)
!					reslt(2,2) = heaviside_w(1,w_index)

!					call fill_fourier_parameters(Trange,wrange,dw,gg)
!					write(*,*) 'Dw,max(Lambda_w)',Delta_w_index*dw, maxval(abs(exc1(1)%dLambda_exc_w(2,1,:)))
!					stop
!				end if

		else if(type == '2') then
			reslt = -(0,1)*(matmul(Hm_2,rho)-matmul(rho,Hm_conjg))
!			do s=1,Nl1
!				reslt = reslt +													&
!					matmul(matmul(exc1(s)%dLambda_2_exc_w(:,:,w_index),rho),		&
!					transpose(conjg(exc1(s)%K_exc)))
!			end do

!			do s=1,Nl2
!				reslt = reslt +													&
!					matmul(matmul(exc2(s)%K_2_exc,rho),							&
!					transpose(conjg(exc2(s)%dLambda_exc_w(:,:,min_w_index))) )
!			end do
		else
			call print_error_message(-1,  'wrong type in fill_superoperator_M_w()')
			stop
		end if

		if(vratit_nulu) then
			reslt = 0.0_dp
		endif

		do i=1, first_N
		do j=1, second_N

		M_w(SUPERINDEX_FROM_K_L(i,j,second_N),SUPERINDEX_FROM_K_L(k,l,second_N)) = &
			reslt(i,j)

		end do
		end do

	end do
	end do


 	end subroutine fill_superoperator_M_w

 	subroutine fill_superoperator_L(LLL,type)
 		character, intent(in) :: type
        complex(dpc), dimension(:,:), intent(out) :: LLL

 		complex(dpc), dimension(N1_from_type(type),N2_from_type(type)) 	:: rho, reslt , smazat
 		complex(dpc), dimension(Nl1,Nl1) 									:: Ham
 		complex(dpc), dimension(Nl2,Nl2) 									:: Ham_2
 		integer(i4b) 			:: i,j,k,l,s, first_N, second_N

 	if(	(.not.(size(LLL,1) == N1_from_type(type)*N2_from_type(type))) .and. &
 		(.not.(size(LLL,2) == N1_from_type(type)*N2_from_type(type))) ) then
 			call print_error_message(-1,  'Dimension error in fill_superoperator_L()')
 			stop
 	end if

	call prepare_Hamiltonian_exc(Ham, 'E');
	if(use_twoexcitons) then
		call prepare_Hamiltonian_exc(Ham_2,'F');
	end if

	first_N 	= N1_from_type(type)
	second_N	= N2_from_type(type)

	do k=1, first_N
	do l=1, second_N
		rho = 0.0_dp
		rho(k,l) = 1.0_dp

		if(type == 'O') then
			reslt = matmul(Ham,rho)
		else if(type == 'E') then
			reslt = matmul(Ham,rho) - matmul(rho,Ham)
		else if(type == '2') then
			reslt = matmul(Ham_2,rho) + matmul(rho,Ham)
		else
			call print_error_message(-1,  'wrong type in fill_superoperator_L()')
			stop
		end if

		do i=1, first_N
		do j=1, second_N

		LLL(SUPERINDEX_FROM_K_L(i,j,second_N),SUPERINDEX_FROM_K_L(k,l,second_N)) = &
			reslt(i,j)

		end do
		end do
	end do
	end do

 	end subroutine fill_superoperator_L

 	subroutine fill_fourier_parameters(Trange,wrange,dw,gg)
 		real(dp), intent(out) :: Trange,wrange,dw,gg

 		Trange 	= dt*grid_Nt*2
		wrange	= 2*PI*grid_Nt*2/Trange
		dw		= wrange/(grid_Nt*2)
		gg 		= 1/NOSE_QME_NZ_REG_TIME_PER_TMAX/(Trange/2)
 	end subroutine


















!------------------------------------- QME global --------------------------------

	! 'Global' functions
 	subroutine init_qme(submethod)
		character, intent(in) 	:: submethod
		integer					:: i,j,k

		if(.not.(											&
			submethod == 'Q' .or. submethod == 'q' .or. 	&
			submethod == 'r' .or. submethod == 's' .or.	&
			submethod == 'R' .or. submethod == 'S' .or.	&
			submethod == 'P' .or. submethod == 'p' .or.	&
			submethod == 'u' .or. submethod == 'U' .or. &
			submethod == '#'			&
			)) then
			call print_error_message(-1,'wrong submethod in init_qme()')
			stop
		end if

		ALLOCATE(heaviside_w, (1,grid_Nt*2) )

		do i=1,grid_Nt
			heaviside_w(1,i) = 1.0_dp
		end do
		do i=grid_Nt+1,grid_Nt*2
			heaviside_w(1,i) = 0.0_dp
		end do

		call fft_row_regularized(heaviside_w, 1, NOSE_QME_NZ_REG_TIME_PER_TMAX)


		Nl1 = iblocks(1,1)%eblock%N1
		Nl2 = iblocks(1,1)%eblock%N1*(iblocks(1,1)%eblock%N1-1)/2

 		allocate(exc1(Nl1))

 		do i=1, Nl1
 			! 'E' block of 1-exciton states
 			ALLOCATE(exc1(i)%K_exc, (Nl1,Nl1))
 			if(submethod == 'Q' .or. submethod == 'q') then
 				ALLOCATE(exc1(i)%dLambda_exc, (Nl1,Nl1,grid_Nt*2))
 				ALLOCATE(exc1(i)%dLambda_exc_w, (Nl1,Nl1,grid_Nt*2))
 			else if(submethod == 'r' .or. submethod == 's' .or. &
 					submethod == 'R' .or. submethod == 'S' .or. &
 					submethod == 'P' .or. submethod == 'p' .or. &
 					submethod == '#') then
 				ALLOCATE(exc1(i)%dLambda_exc, (Nl1,Nl1,grid_Nt))
 				ALLOCATE(exc1(i)%Lambda_exc, (Nl1,Nl1,grid_Nt))
 			else if(submethod == 'u' .or. submethod == 'U') then
 				ALLOCATE(exc1(i)%dLambda_exc, (Nl1,Nl1,grid_Nt))
 				ALLOCATE(exc1(i)%Lambda_exc, (Nl1,Nl1,grid_Nt))
 			end if

 		end do

 		if(use_twoexcitons) then
 			allocate(exc2(Nl2))

 			do i=1, Nl2

 				! 'F' block of 2-exciton states
 				ALLOCATE(exc2(i)%K_2_exc, (Nl2,Nl2))

 				if(submethod == 'Q' .or. submethod == 'q') then
 					ALLOCATE(exc2(i)%dLambda_2_exc, (Nl2,Nl2,grid_Nt*2))
 					ALLOCATE(exc2(i)%dLambda_2_exc_w, (Nl2,Nl2,grid_Nt*2))
 				else if(submethod == 'r' .or. submethod == 's' .or. &
 					 	submethod == 'R' .or. submethod == 'S' .or. &
 					 	submethod == 'P' .or. submethod == 'p' .or. &
 					 	submethod == '#') then
 					ALLOCATE(exc2(i)%dLambda_2_exc, (Nl2,Nl2,grid_Nt))
 					ALLOCATE(exc2(i)%Lambda_2_exc, (Nl2,Nl2,grid_Nt))
 				else if(submethod == 'u' .or. submethod == 'U') then
 					ALLOCATE(exc2(i)%dLambda_2_exc, (Nl2,Nl2,grid_Nt))
 					ALLOCATE(exc2(i)%Lambda_2_exc, (Nl2,Nl2,grid_Nt))
 				end if
	 		end do

 		end if

 	end subroutine init_qme

 	subroutine deinit_qme(submethod)
 		character, intent(in)	:: submethod
		integer					:: i,j,k

		if(.not.(											&
			submethod == 'Q' .or. submethod == 'q' .or. 	&
			submethod == 'r' .or. submethod == 's' .or.	&
			submethod == 'R' .or. submethod == 'S' .or.	&
			submethod == 'P' .or. submethod == 'p' .or.	&
			submethod == 'u' .or. submethod == 'U'	.or. &
			submethod == '#'		&
			)) then
			call print_error_message(-1,'wrong submethod in deinit_qme()')
			stop
		end if
		DEALLOCATE(heaviside_w)

 		do i=1, Nl1
 			! 'E' block of 1-exciton states
 			DEALLOCATE(exc1(i)%K_exc)
 			DEALLOCATE(exc1(i)%dLambda_exc)
 			if(submethod == 'Q' .or. submethod == 'q') then
 				DEALLOCATE(exc1(i)%dLambda_exc_w)
 			else if(submethod == 'r' .or. submethod == 's' .or. &
 					submethod == 'R' .or. submethod == 'S' .or. &
					submethod == 'P' .or. submethod == 'p' .or. &
					submethod == '#') then
 				DEALLOCATE(exc1(i)%Lambda_exc)
 			else if(submethod == 'u' .or. submethod == 'U') then
 				DEALLOCATE(exc1(i)%Lambda_exc)
 			end if
 		end do

 		if(use_twoexcitons) then
 			do i=1, Nl2
				! 'F' block of 2-exciton states
	 			DEALLOCATE(exc2(i)%K_2_exc)
 				DEALLOCATE(exc2(i)%dLambda_2_exc)
 				if(submethod == 'Q' .or. submethod == 'q') then
 					DEALLOCATE(exc2(i)%dLambda_2_exc_w)
 				else if(submethod == 'r' .or. submethod == 's' .or. &
 						submethod == 'R' .or. submethod == 'S' .or. &
 						submethod == 'P' .or. submethod == 'p' .or. &
 						submethod == '#') then
 					DEALLOCATE(exc2(i)%Lambda_2_exc)
 				else if(submethod == 'u' .or. submethod == 'U') then
 					DEALLOCATE(exc2(i)%Lambda_2_exc)
 				end if
 			end do

 			deallocate(exc2)
 		end if

 		 deallocate(exc1)
 	end subroutine deinit_qme

	subroutine perform_qme(rho, rho0, type, submethod)
		character, intent(in) 		:: submethod, type
		complex(dpc), dimension(:,:), intent(in)	:: rho0
		complex(dpc), dimension(:,:,:), intent(out)	:: rho

		complex(dpc), dimension(	N1_from_type(type)*N2_from_type(type),		&
 									N1_from_type(type)*N2_from_type(type) ) 	&
 					   				:: superoperator_M_w, superoperator_L, 		&
 					   				   superoperator_G, superoperator_1, superoperator_G_inv
 		integer(i4b) 				:: i,j,k, w_index

 		real(dp)					:: gg, Trange, wrange, dw, omega, time

		if(N1_from_type(type) < 1 .or. N2_from_type(type) < 1) then
			call print_error_message(-1,  'N1_from_type or N2_from_type given wrong type in perform_qme()')
			stop
		end if

		if(.not.(N1_from_type(type) == size(rho0,1) .and. N2_from_type(type) == size(rho0,2))) then
			call print_error_message(-1,  'rho0 of wrong dimension in perform_qme()')
			stop
		end if

		if(.not.(N1_from_type(type) == size(rho,1) .and. N2_from_type(type) == size(rho,2) &
				.and. size(rho,3) == grid_Nt*2 )) then
			call print_error_message(-1,  'rho of wrong dimension in perform_qme()')
			stop
		end if

		if(.not.(submethod == 'q' .or. submethod == 'Q')) then
			call print_error_message(-1, 'incorrect submethod in perform_qme()')
			stop
		end if


		superoperator_1 = 0.0_dp

		do i=1, size(superoperator_1,1)
			superoperator_1(i,i) = 1.0_dp
		end do

		call fill_fourier_parameters(Trange,wrange,dw,gg)

		call prepare_K_excs()
		call prepare_dLambdas()
		call fourier_dLambdas()
		call generate_warnings(submethod)

		call fill_superoperator_L(superoperator_L, type)


		do w_index=1, size(exc1(1)%dLambda_exc_w,3)
			if(w_index <= grid_Nt*2/2) then
				omega = (w_index-1)*dw
			else
				omega = ((w_index-1)-grid_Nt*2)*dw
			end if

			call fill_superoperator_M_w(superoperator_M_w, w_index, type)

			if(submethod == 'Q') then
				call secularize_superoperator(superoperator_M_w, type)
			end if

			superoperator_G_inv = superoperator_1*(-omega*(0.0_dp,1.0_dp) + gg) &
								- dt*superoperator_M_w + (0.0_dp,1.0_dp)*superoperator_L

			call inv(superoperator_G_inv, superoperator_G)

			!call perform_superoperator(superoperator_M_w*dt,rho0/dt,rho(:,:,w_index))
			call perform_superoperator(superoperator_G,rho0/dt,rho(:,:,w_index))


		end do

		call fft_row_regularized(rho, -1, NOSE_QME_NZ_REG_TIME_PER_TMAX)


	end subroutine perform_qme




















































!------------------------------------- REDFIELD --------------------------------

 	subroutine integrate_Lambdas()
 		integer(i4b) 		:: i,j,k

		! 'E' block of 1-exciton states
 		do i=1, Nl1
			exc1(i)%Lambda_exc(:,:,1) = 0.0_dp
			do k=1, grid_Nt-1
				exc1(i)%Lambda_exc(:,:,k+1) = exc1(i)%Lambda_exc(:,:,k)	+				&
					dt*(exc1(i)%dLambda_exc(:,:,k)+exc1(i)%dLambda_exc(:,:,k+1))/2.0_dp
			end do
		end do

		if(use_twoexcitons) then
			! 'F' block of 2-exciton states
			do i=1, Nl2
				exc2(i)%Lambda_2_exc(:,:,1) = 0.0_dp
				do k=1, grid_Nt-1
					exc2(i)%Lambda_2_exc(:,:,k+1) = exc2(i)%Lambda_2_exc(:,:,k)	+			&
						dt*(exc2(i)%dLambda_2_exc(:,:,k)+exc2(i)%dLambda_2_exc(:,:,k+1))/2.0_dp
				end do
			end do
		end if

 	end subroutine integrate_Lambdas


 	subroutine fill_superoperator_R(RRR,t_index,type)
 		character, intent(in) :: type
        integer, intent(in) :: t_index
        complex(dpc), dimension(:,:), intent(out) :: RRR

 		complex(dpc), dimension(N1_from_type(type),N2_from_type(type)) 	&
 											:: rho, reslt
 		complex(dpc), dimension(Nl1,Nl1)	:: Hm, Hm_conjg
 		complex(dpc), dimension(Nl2,Nl2)	:: Hm_2, Hm_2_conjg
 		integer(i4b) 	:: i,j,k,l,s, first_N, second_N

 	if(	(.not.(size(RRR,1) == N1_from_type(type)*N2_from_type(type))) .and. &
 		(.not.(size(RRR,2) == N1_from_type(type)*N2_from_type(type))) ) then
 			call print_error_message(-1,  'Dimension error in fill_superoperator_R()')
 			stop
 	end if

 	if(	(t_index > size(exc1(1)%Lambda_exc,3) ) .or. &
 		(t_index < 1) ) then
 			call print_error_message(-1,  'Time/frequency dimension error in fill_superoperator_R()')
 			stop
 	end if



	first_N 	= N1_from_type(type)
	second_N	= N2_from_type(type)

	do k=1, first_N
	do l=1, second_N
		rho = 0.0_dp
		rho(k,l) = 1.0_dp


		Hm 			= 0.0_dp
		Hm_2 		= 0.0_dp

		Hm_conjg	= 0.0_dp
		Hm_2_conjg	= 0.0_dp

		! K_exc from 'F' and K_2_exc from 'E' are zero, only terms from
		! particular block contribute

		do i = 1, Nl1

		! 'E' block of 1-exciton states
			Hm 		= Hm 		&
		+ (0,-1)*matmul(exc1(i)%K_exc,exc1(i)%Lambda_exc(:,:,t_index))

			Hm_conjg= Hm_conjg 	&
		+ transpose(conjg( 		&
		  (0,-1)*matmul(exc1(i)%K_exc,exc1(i)%Lambda_exc(:,:,t_index)) &
		  ))

		end do


		if(type == '2') then

		! 'F' block of 1-exciton states		--- no contrinution

		do i = 1, Nl2

		! 'E' block of 2-exciton states     --- no contribution

		! 'F' block of 2-exciton states
			Hm_2 		= Hm_2			&
		+ (0,-1)*matmul(exc2(i)%K_2_exc,exc2(i)%Lambda_2_exc(:,:,t_index))
			Hm_2_conjg 	= Hm_2_conjg	&
		+ transpose(conjg( 				&
		  (0,-1)*matmul(exc2(i)%K_2_exc,exc2(i)%Lambda_2_exc(:,:,t_index)) &
		  ))

		end do
		end if






		if(type == 'O') then
			reslt = -(0,1)*(matmul(Hm,rho) )
		else if(type == 'E') then
			reslt = -(0,1)*(matmul(Hm,rho) - matmul(rho,Hm_conjg))
!			reslt = 0.0

			do s=1,Nl1
				reslt = reslt 													&
					+ matmul(matmul(exc1(s)%K_exc,rho),							&
					transpose(conjg(exc1(s)%Lambda_exc(:,:,t_index))) ) 		&
					+ matmul(matmul(exc1(s)%Lambda_exc(:,:,t_index),rho),		&
					transpose(conjg(exc1(s)%K_exc)) )

			end do

		else if(type == '2') then
			reslt = -(0,1)*(matmul(Hm_2,rho)-matmul(rho,Hm_conjg))
!			do s=1,Nl1
!				reslt = reslt +													&
!					matmul(matmul(exc1(s)%Lambda_2_exc(:,:,t_index),rho),		&
!					transpose(conjg(exc1(s)%K_exc)))
!			end do

!			do s=1,Nl2
!				reslt = reslt +													&
!					matmul(matmul(exc2(s)%K_2_exc,rho),							&
!					transpose(conjg(exc2(s)%Lambda_exc(:,:,t_index))) )
!			end do
		else
			call print_error_message(-1,  'wrong type in fill_superoperator_R()')
			stop
		end if

		do i=1, first_N
		do j=1, second_N

		RRR(SUPERINDEX_FROM_K_L(i,j,second_N),SUPERINDEX_FROM_K_L(k,l,second_N)) = &
			reslt(i,j)

		end do
		end do

	end do
	end do


 	end subroutine fill_superoperator_R

 	subroutine calculate_derivative_redf(t, y, dydt)
 		real(dp), intent(in) :: t
 		complex(dpc), dimension(:,:), intent(in) 	:: y
 		complex(dpc), dimension(:,:), intent(out) 	:: dydt

 		integer(i4b) :: t_index, index_i, index_j

		complex(dpc), dimension(	N1_from_type(global_type)*N2_from_type(global_type),	&
									N1_from_type(global_type)*N2_from_type(global_type) ) 	&
				   				:: superoperator_R
		complex(dpc), dimension(	N1_from_type(global_type)*N2_from_type(global_type),	&
									N1_from_type(global_type)*N2_from_type(global_type) ) 	&
				   				:: whole_superoperator
		complex(dpc), dimension(	N1_from_type(global_type) )   &
				   				:: eigenvalues
		complex(dpc), dimension(	N1_from_type(global_type),N1_from_type(global_type) )   &
				   				:: populations_only_superoperator, eigenvectors, eigenvalueMatrix

		character(len=256)		:: buff

		if(	(.not. size(y,1) == N1_from_type(global_type)).or. &
			(.not. size(y,2) == N2_from_type(global_type)).or. &
			(.not. size(dydt,1) == N1_from_type(global_type)).or. &
			(.not. size(dydt,2) == N2_from_type(global_type))	) then
			call print_error_message(-1, "dimension error in calculate_derivative_redf")
			stop
		end if

		t_index = INT(t/dt)+1

		if(t_index < 1 .or. t_index > grid_Nt) then
			call print_error_message(-1, "time out of range in calculate_derivative_redf")
		end if

		if(global_submethod == 'r' .or. global_submethod == 's' .or. &
		   global_submethod == 'p') then
			call fill_superoperator_R(superoperator_R, t_index, global_type)
		elseif(global_submethod == 'R' .or. global_submethod == 'S' .or. global_submethod == 'P') then
			call fill_superoperator_R(superoperator_R, grid_Nt, global_type)
		elseif(global_submethod == '#') then
			superoperator_R = 0.0_dp
		endif

		if((global_submethod == 'p' .or. global_submethod == 'P') .and. global_type == 'E') then
			call phenom_relax(superoperator_R)
		end if

		if(global_submethod == 's' .or. global_submethod == 'S' .or. &
		   global_submethod == 'p' .or. global_submethod == 'P') then
			call secularize_superoperator(superoperator_R, global_type)
		end if

		! This is quite minor part. All elements of whole superoperator that do
		! not directly affect population space are removed and we calculate the
		! eigenvalues of the rest. If the result has positive eigenvalues, not
		! only the positivity is lost, but solution can diverge. This occurs
		! even for secular approximation, though not for long times. This part
		! basicaly writes warning to be careful in these cases - they should't be
		! trusted even in secular approximation.
		if(t_index == grid_Nt-1 .and. global_type == 'E') then
			populations_only_superoperator = 0
			eigenvectors = 0
			eigenvalueMatrix = 0

			whole_superoperator = ((0.0,-1.0)*global_superop_L+superoperator_R)

			do index_i=1,N1_from_type(global_type)**2,N1_from_type(global_type)+1
			do index_j=1,N1_from_type(global_type)**2,N1_from_type(global_type)+1

			populations_only_superoperator((index_i-1)/(N1_from_type(global_type)+1)+1,&
				 (index_j-1)/(N1_from_type(global_type)+1)+1) = whole_superoperator(index_i,index_j)

			end do
			end do

			call spec(populations_only_superoperator,eigenvectors,eigenvalueMatrix)

			do index_i=1,N1_from_type(global_type)
				eigenvalues(index_i) = eigenvalueMatrix(index_i,index_i)
			end do

			if( maxval(real(eigenvalues)) > 1.0e-5 ) then
				call print_warning_message('Redfield tensor has positive population eigenvalues at big times',3)
				write(buff,*) 'Re(eig)=',real(eigenvalues)
				call print_warning_message(trim(buff),3)
			end if
		end if
		! end of eigenvalues writing part

		call perform_superoperator((0.0,-1.0)*global_superop_L+superoperator_R, y, dydt)

 	end subroutine calculate_derivative_redf

	subroutine phenom_relax(redfieldsop)

		complex(dpc), dimension(:,:), intent(inout) :: redfieldsop
		complex(dpc), dimension(size(redfieldsop,1),size(redfieldsop,1)) :: RelaxGamma ! relax. superoperator
		integer :: i

		redfieldsop = redfieldsop - phenom_Rates 

	end subroutine phenom_relax


	subroutine prepare_phenom_Rates()

		real(dp) :: gammaA, gammaB ! temporary means: pop. relaxation elements set by hand
		real(dp) :: kappa ! for testing
		gammaA = 1.0_dp/2000000_dp		! e.g. chlorophyl's excitation life-time ~ 2 ns
		gammaB = 1.0_dp/5000_dp			! e.g. carotenoid's excitation life-time ~ 5ps

		! relaxation superoperator in site representation
		phenom_Rates = 0.0_dp
		phenom_Rates(1,1) = gammaA; phenom_Rates(4,4) = gammaB ! population relaxation
		phenom_Rates(2,2) = (0.5_dp)*(gammaA + gammaB) ! coherence dephasing
		phenom_Rates(3,3) = phenom_Rates(2,2)

		call superops_to_exc(phenom_Rates, 'E')

		open(unit=11,file=trim(file_join(out_dir,"phenom_rates_exc.dat")))
			write(11,'(4f16.8)') real(phenom_Rates(1,1)), real(phenom_Rates(1,2)), real(phenom_Rates(1,3)), real(phenom_Rates(1,4))
			write(11,'(4f16.8)') real(phenom_Rates(2,1)), real(phenom_Rates(2,2)), real(phenom_Rates(2,3)), real(phenom_Rates(2,4))
			write(11,'(4f16.8)') real(phenom_Rates(3,1)), real(phenom_Rates(3,2)), real(phenom_Rates(3,3)), real(phenom_Rates(3,4))
			write(11,'(4f16.8)') real(phenom_Rates(4,1)), real(phenom_Rates(4,2)), real(phenom_Rates(4,3)), real(phenom_Rates(4,4))
		close(unit=11)


	end subroutine prepare_phenom_Rates


	subroutine perform_redfield(rho, rho0, type, submethod)

		character, intent(in) 		:: submethod, type
		complex(dpc), dimension(:,:), intent(in)	:: rho0
		complex(dpc), dimension(:,:,:), intent(out)	:: rho

		complex(dpc), dimension(size(rho,1),size(rho,2))	:: drho
		real(dp) :: t, tt
		integer(i4b) :: i, t_index


		if(N1_from_type(type) < 1 .or. N2_from_type(type) < 1) then
			call print_error_message(-1,  'N1_from_type or N2_from_type given wrong type in perform_redfield()')
			stop
		end if

		if(.not.(N1_from_type(type) == size(rho0,1) .and. N2_from_type(type) == size(rho0,2))) then
			call print_error_message(-1,  'rho0 of wrong dimension in perform_redfield()')
			stop
		end if

		if(.not.(N1_from_type(type) == size(rho,1) .and. N2_from_type(type) == size(rho,2) &
				.and. size(rho,3) >= grid_Nt )) then
			call print_error_message(-1,  'rho of wrong dimension in perform_redfield()')
			stop
		end if

		if(.not.(submethod == 'r' .or. submethod == 's' .or. &
				  submethod == 'R' .or. submethod == 'S' .or. &
				  submethod == 'P' .or. submethod == 'p' .or. &
				  submethod == 'u' .or. submethod == 'U' .or. &
				  submethod == '#') ) then
			call print_error_message(-1, 'incorrect submethod in perform_redfield()')
			stop
		end if

		call prepare_K_excs()
		call prepare_dLambdas()
		call integrate_Lambdas()

		if(submethod == 'p' .or. submethod == 'P') then
			ALLOCATE(phenom_Rates,(N1_from_type('E')*N1_from_type('E'),N1_from_type('E')*N1_from_type('E')))
			call prepare_phenom_Rates
		end if

		if(submethod == 'u' .or. submethod == 'U') then
!			call prepare_Upsilonas()
		end if
		call generate_warnings(submethod)



		global_type = type
		global_submethod = submethod
		ALLOCATE(global_superop_L,(N1_from_type(type)*N2_from_type(type),N1_from_type(type)*N2_from_type(type)))
		call fill_superoperator_L(global_superop_L, global_type)

		if(global_type == 'O' .and. .not. Redfield_printed) then
			call write_redfield_tensor2('O')
			Redfield_printed = .true.
			write(*,*) 'Redfield printed'
		end if

		rho(:,:,1) = rho0
		drho = 0.0_dp

		do t_index = 1, grid_Nt-1
			tt = (t_index-1)*dt

			if(submethod == 'r' .or. submethod == 's' .or. &
			   submethod == 'P' .or. submethod == 'p' .or. &
			   submethod == 'R' .or. submethod == 'S' .or. &
			   submethod == '#') then
				call calculate_derivative_redf(tt,rho(:,:,t_index),drho)
				call ode_rk4(rho(:,:,t_index),drho,tt,dt,rho(:,:,t_index+1),calculate_derivative_redf)

			else if(submethod == 'u' .or. submethod == 'U') then
				call calculate_derivative_red_tau(tt,rho(:,:,t_index),drho)
				call ode_rk4(rho(:,:,t_index),drho,tt,dt,rho(:,:,t_index+1),calculate_derivative_red_tau)

			end if

		end do

		DEALLOCATE(global_superop_L)
		
		if(allocated(phenom_Rates)) then
			DEALLOCATE(phenom_Rates)
		end if

		if(global_type == 'E' .and. submethod == '#') then
			! only 12-excitonic coherence of a dimer is calculated
			do t_index = 2, grid_Nt
				tt = (t_index-1)*dt
				if(tt+tau_of_projector > (grid_Nt-1)*dt) then
					exit
				end if

!				rho(1,2,t_index) = &!rho(1,2,t_index)* &
!				exp( -goft_exciton(2,2,tt+tau_of_projector)-conjg(goft_exciton(1,1,tt))  )   !/&
!				!exp( -goft_exciton(2,2,0.0+tau_of_projector)-conjg(goft_exciton(1,1,0.0_dp)) )

				if(abs(rho(1,2,2)) > 0) then
					rho(1,2,t_index) = rho(1,2,t_index)* &
					exp(+goft_exciton(2,2,tau_of_projector))/(rho(1,2,2)/abs(rho(1,2,2)) )* &
					exp( -conjg(goft_exciton(2,2,tt+tau_of_projector))  -goft_exciton(1,1,tt)  )
				end if

			end do

		end if

	end subroutine perform_redfield

























!------------------------------------- tau-Projectors REDFIELD --------------------------------


	function goft_exciton(m,n,t) result(goft)
		real(dp), intent(in)		:: t
		integer(i4b), intent(in)	:: m,n
		complex(dpc)				:: goft
		integer						:: a, t_index

 		t_index = INT(t/dt)+1

		goft = 0.0_dp

		do a=1,Nl1
	 		if(t_index < 1 .or. t_index > size(all_goft(iblocks(1,1)%sblock%gindex(a))%gg,1)) then
 				call print_error_message(-1,"tau exceeds goft size in goft_exciton()")
 				stop
 			end if

			goft = goft +		 													&
        		(iblocks(1,1)%eblock%SS(a,m)**2)*(iblocks(1,1)%eblock%SS(a,n)**2) 	&
				*all_goft(iblocks(1,1)%sblock%gindex(a))%gg(t_index)
		end do

	end function goft_exciton

	function dgoft_exciton(m,n,t) result(dgoft)
		real(dp), intent(in)		:: t
		integer(i4b), intent(in)	:: m,n
		complex(dpc)				:: dgoft
		integer						:: a, t_index

 		t_index = INT(t/dt)+1

		dgoft = 0.0_dp

		do a=1,Nl1
 			if(t_index < 1 .or. t_index > size(all_hoft(iblocks(1,1)%sblock%gindex(a))%gg,1)) then
 				call print_error_message(-1,"tau exceeds goft size in dgoft_exciton()")
 				stop
 			end if

			dgoft = dgoft +			 												&
        		(iblocks(1,1)%eblock%SS(a,m)**2)*(iblocks(1,1)%eblock%SS(a,n)**2) 	&
				*all_hoft(iblocks(1,1)%sblock%gindex(a))%gg(t_index)
		end do

	end function dgoft_exciton

	function dgoft_half_exciton(m,n,t) result(dgoft)
		real(dp), intent(in)		:: t
		integer(i4b), intent(in)	:: m,n
		complex(dpc)				:: dgoft
		integer						:: t_index

 		t_index = INT(t/dt)+1

		dgoft = 0.0_dp

		if(t_index < 1 .or. t_index > size(all_hoft(iblocks(1,1)%sblock%gindex(m))%gg,1)) then
			call print_error_message(-1,"tau exceeds goft size in dgoft_half_exciton()")
			stop
		end if

		dgoft = dgoft +			 												&
        	(iblocks(1,1)%eblock%SS(m,n)**2) 									&
			*all_hoft(iblocks(1,1)%sblock%gindex(m))%gg(t_index)


	end function dgoft_half_exciton

	function dgoft_site(m,n,t) result(dgoft)
		real(dp), intent(in)		:: t
		integer(i4b), intent(in)	:: m,n
		complex(dpc)				:: dgoft
		integer	(i4b)				:: a, t_index

 		t_index = INT(t/dt)+1

		dgoft = 0.0_dp

		if(m == n) then

			if(m < 1 .or. m > size(iblocks(1,1)%sblock%gindex)) then
 				call print_error_message(-1,"m exceeds Nl1 size in dgoft_site()")
 				stop
 			end if

			if(t_index < 1 .or. t_index > size(all_hoft(iblocks(1,1)%sblock%gindex(m))%gg,1)) then
 				call print_error_message(-1,"tau exceeds goft size in dgoft_site()")
 				stop
 			end if

			dgoft = all_hoft(iblocks(1,1)%sblock%gindex(m))%gg(t_index)

		end if

	end function dgoft_site

! 	subroutine weigth_superoperator_by_exc_goft(superop, t)
 !		complex(dpc), dimension(:,:), intent(inout)	:: superop
 !		real(dp), intent(in)							:: t
!
 !		integer			:: i,j,k,l, tau_index
 !		complex(dpc), dimension(size(superop,1),size(superop,2))		&
 !				:: Ue_superop_exc, Ue1_superop_exc
!
 !		tau_index = INT(tau/dt)+1
!
 !		if(tau_index < 1 .or. tau_index > size(iblocks(1,1)%eblock%gg,2)) then
 !			call print_error_message(-1,"tau exceeds goft size in weigth_superoperator_by_exc_goft()")
 !			stop
 !		end if
!
 !		if(tau_index < 1 .or. tau_index > size(iblocks(1,1)%eblock%gg,2)) then
 !			call print_error_message(-1,"tau exceeds goft size in weigth_superoperator_by_exc_goft()")
 !			stop
 !		end if
!
!		call fill_superoperator_L(Ue_superop_exc, 'E')
!
!		do i=1,size(Ue_superop_exc,1)
!			Ue_superop_exc(i,i) = exp(- (0.0,1.0) * Ue_superop_exc(i,i))
!		end do
!
!		Ue1_superop_exc = conjg(Ue_superop_exc) ! transposition ommited for diagonal matrix
!
!		superop = matmul(matmul(Ue1_superop_exc,superop),Ue_superop_exc)
!
 !		do i = 1, size(superop, 1)
 !		do j = 1, size(superop, 2)
 !			superop(i,j) = superop(i,j) * 																	&
!				exp(-conjg(goft_exciton(L_FROM_SUPERINDEX(j,Nl1),L_FROM_SUPERINDEX(j,Nl1),tau) )) * 	&
!				exp(conjg(goft_exciton(L_FROM_SUPERINDEX(i,Nl1),L_FROM_SUPERINDEX(i,Nl1),tau) ))
 !		end do
 !		end do
!
 !		superop = matmul(matmul(Ue_superop_exc,superop),Ue1_superop_exc)
!
! 	end subroutine weigth_superoperator_by_exc_goft

 	subroutine fill_superoperator_R_tau_addition(Rta,t_index)
        integer, intent(in) :: t_index
        complex(dpc), dimension(:,:), intent(out) :: Rta

  		complex(dpc), dimension(Nl1,Nl1) 	&
 											:: rho, reslt
 		complex(dpc), dimension(Nl1,Nl1)	:: Ue_exc, Ue1_exc, K_tilde, K_tilde_tni
 		integer(i4b)	:: i,j,m,l,u,v
 		real(dp)		:: t

 	if(	(.not.(size(Rta,1) == Nl1*Nl1)) .and. &
 		(.not.(size(Rta,2) == Nl1*Nl1)) ) then
 			call print_error_message(-1,  'Dimension error in fill_superoperator_R_tau_addition()')
 			stop
 	end if

 	if(	(t_index > size(exc1(1)%Lambda_exc,3) ) .or. &
 		(t_index < 1) ) then
 			call print_error_message(-1,  'Time dimension error in fill_superoperator_R_tau_addition()')
 			stop
 	end if

	t = (t_index-1)*dt

	! if we are out of covered goft data, we assume that the value is the
	! same as for t_index = max_t_index
	if(.not. ( t+tau_of_projector <= dt*(grid_Nt-1) )) then
		t = (grid_Nt-1)*dt - tau_of_projector
	end if

	call prepare_Ue_exc(Ue_exc, t, 'E')
	Ue1_exc = conjg(transpose(Ue_exc))


	do u=1, Nl1
	do v=1, Nl1
		rho = 0.0_dp
		rho(u,v) = 1.0_dp

		reslt = 0.0

		do m=1,Nl1
		do l=1,Nl1

!		HALF-EXCITON
!			K_tilde 			= 0.0_dp
!			K_tilde(m,m)	 	= 1.0_dp
			K_tilde 			= exc1(m)%K_exc

			K_tilde_tni 		= 0.0_dp
			K_tilde_tni(l,l)	= 1.0_dp
!			K_tilde_tni 		= exc1(l)%K_exc
			K_tilde_tni 		= matmul(matmul(Ue_exc,K_tilde_tni),Ue1_exc)



			reslt = reslt + (matmul(matmul(K_tilde,rho),K_tilde_tni)  	&
					      - matmul(matmul(rho,K_tilde_tni),K_tilde) ) *	&
					conjg(dgoft_half_exciton(m,l,t+tau_of_projector)-dgoft_half_exciton(m,l,t) )


!		SITE - doesn't work - and is not expected to work
!			K_tilde 			= exc1(m)%K_exc
!
!			K_tilde_tni 		= exc1(m)%K_exc
!			K_tilde_tni 		= matmul(matmul(Ue_exc,K_tilde_tni),Ue1_exc)
!
!			reslt = reslt + ( matmul(K_tilde,matmul(rho,K_tilde_tni))  	&
!					      - matmul(matmul(rho,K_tilde_tni),K_tilde) ) *	&
!					conjg(dgoft_site(m,l,t+tau_of_projector)-dgoft_site(m,l,t))

!			reslt = reslt * exp(conjg(goft_exciton(l,l,tau_of_projector) ))

		end do
		end do

		do i=1, Nl1
		do j=1, Nl1

!		reslt(i,j) = reslt(i,j) * exp(-conjg(goft_exciton(j,j,tau_of_projector) ))

		Rta(SUPERINDEX_FROM_K_L(i,j,Nl1),SUPERINDEX_FROM_K_L(u,v,Nl1)) = &
			reslt(i,j)

		end do
		end do

	end do
	end do



 	end subroutine fill_superoperator_R_tau_addition

 	subroutine calculate_derivative_red_tau(t, y, dydt)
 		real(dp), intent(in) :: t
 		complex(dpc), dimension(:,:), intent(in) 	:: y
 		complex(dpc), dimension(:,:), intent(out) 	:: dydt

 		integer(i4b) :: t_index

		complex(dpc), dimension(	N1_from_type(global_type)*N2_from_type(global_type),	&
									N1_from_type(global_type)*N2_from_type(global_type) ) 	&
				   				:: superoperator_R, superoperator_Rta

		if(	(.not. size(y,1) == N1_from_type(global_type)).or. &
			(.not. size(y,2) == N2_from_type(global_type)).or. &
			(.not. size(dydt,1) == N1_from_type(global_type)).or. &
			(.not. size(dydt,2) == N2_from_type(global_type))	) then
			call print_error_message(-1, "dimension error in calculate_derivative_red_tau")
			stop
		end if

		t_index = INT(t/dt)+1

		if(t_index < 1 .or. t_index > grid_Nt) then
			call print_error_message(-1, "time out of range in calculate_derivative_red_tau")
		end if

		call fill_superoperator_R(superoperator_R, t_index, global_type)
		!call weigth_superoperator_by_exc_goft(superoperator_R, t)

		call fill_superoperator_R_tau_addition(superoperator_Rta, t_index)

		if(global_submethod == 'U') then
			call secularize_superoperator(superoperator_R, global_type)
			call secularize_superoperator(superoperator_Rta, global_type)
		end if

!		write(*,*) maxval(abs(superoperator_R)), maxval(abs(global_superop_L))

		call perform_superoperator((0.0,-1.0)*global_superop_L+superoperator_R+superoperator_Rta, y, dydt)

	end subroutine calculate_derivative_red_tau

































!------------------------------------- Lambda_w-calc --------------------------------


!	subroutine fill_c_w()
!		double precision :: t
!		integer(i4b) :: i,j,k, N_use, ind1, ind2, ind
!
!		if(type == 'E') then
!			N_use = Nl1
!		else if(type == 'F') then
!			N_use = Nl2
!			ind1 = iblocks(1,1)%eblock%ione1(ind)
!			ind2 = iblocks(1,1)%eblock%ione2(ind)
!		end if
!
!		c_w = 0.0_dp
!		c_2_w = 0.0_dp
!
!		! Lambdas are considered to be zero for t < 0 => maximal index grid_Nt
!		do ind=1,Nl1
!		do i=1,grid_Nt
!			t = (i-1)*dt
!
!			c_w(ind, i) = igofts(iblocks(1,1)%sblock%gindex(ind))%goft%ct(i)
!		end do
!		end do
!
!		do ind=1,Nl2
!		do i=1,grid_Nt
!			t = (i-1)*dt
!
!			c_2_w(ind, i) = (igofts(iblocks(1,1)%sblock%gindex(ind1))%goft%ct(i) + igofts(iblocks(1,1)%sblock%gindex(ind2))%goft%ct(i))
!		end do
!		end do
!
!		end subroutine fill_c_w


























!------------------------------------- PUBLIC --------------------------------

	subroutine fill_evolution_superoperator_nakajima_zwanzig(type, submethod)
		character, intent(in)	:: submethod, type

		complex(dpc), dimension(N1_from_type(type),N2_from_type(type)) 	:: rho0
		complex(dpc), dimension(N1_from_type(type),N2_from_type(type),grid_Nt*2) 	:: rho

		character(len=256) :: cbuff

		integer :: i, k, l, r, s

		logical :: calculate_secular



		if(.not.(type == 'O' .or. type == 'E' .or. type == '2')) then
			call print_error_message(-1, 'error in fill_evolution_superoperator(), wrong type')
			stop
		end if

		if(.not.(	submethod == 'q' .or. submethod == 'Q' .or. &
					submethod == 's' .or. submethod == 'r' .or. &
					submethod == 'S' .or. submethod == 'R' .or. &
					submethod == 'P' .or. submethod == 'p' .or. &
					submethod == 'u' .or. submethod == 'U' .or. &
					submethod == '#')) then
			call print_error_message(-1, 'error in fill_evolution_superoperator(), wrong submethod')
			stop
		end if

		call init_qme(submethod)

		do r=1,N1_from_type(type)
		do s=1,N2_from_type(type)

			calculate_secular = .true.

			! for secular approximation and coherence block, we can run calculation just once, because elements do not mix with each other
			if(submethod == 's' .or. submethod == 'S') then
				calculate_secular = .false.

				if(type == 'O' .or. type == '2') then
					if(r == 1 .and. s == 1) then
						rho 	= 0.0_dp
						calculate_secular = .true.
					end if

					rho0 		= 1.0_dp
				elseif(type == 'E') then
					if(r == s) then
						calculate_secular = .true.
						rho 		= 0.0_dp
						rho0 		= 1.0_dp
!						rho0(r,s) 	= 1.0_dp

						do i=1,N1_from_type(type)
							if(i == s) then
								cycle
							end if

							rho0(i,i) 	= 0.0_dp
						end do
					end if
				end if
			else
				rho0 		= 0.0_dp
				rho0(r,s) 	= 1.0_dp

				rho 		= 0.0_dp
			end if

			if(submethod == 'Q' .or. submethod == 'q') then
				call perform_qme(rho,rho0,type,submethod)
			elseif(submethod == 'r' .or. (submethod == 's' .and. calculate_secular ) .or. &
					submethod == 'R' .or. (submethod == 'S' .and. calculate_secular ) .or. &
					submethod == 'P' .or. submethod == 'p' .or. &
					submethod == 'U' .or. submethod == 'u' .or. &
					submethod == '#') then
				call perform_redfield(rho,rho0,type,submethod)
			end if

			do k=1,N1_from_type(type)
			do l=1,N2_from_type(type)

			rho(k,l,1) = rho0(k,l)

			do i=1,Nt(1)*gt(1),gt(1)
				if(i > Nt(1)*gt(1) .or. i < 1) then
					call print_error_message(-1, 'index error in fill_evolution_superoperator_nakajima_zwangig')
				end if

				if(type == 'O') then
					evops(1,1)%Ueg(k,l, r,s, INT((i-1)/gt(1)) + 1) = rho(k,l,i)

					if((submethod == 's' .or. submethod == 'S') .and. (k /= r .or. l /= s)) then
						evops(1,1)%Ueg(k,l, r,s, INT((i-1)/gt(1)) + 1) = 0.0_dp
					end if
				elseif(type == 'E') then
					evops(1,1)%Uee(k,l, r,s, INT((i-1)/gt(1)) + 1) = rho(k,l,i)

					if((submethod == 's' .or. submethod == 'S') .and. (((k /= r .or. l /= s) .and. k /= l) .or. (k == l .and. r /= s)) ) then
						evops(1,1)%Uee(k,l, r,s, INT((i-1)/gt(1)) + 1) = 0.0_dp
					end if
				elseif(type == '2') then
					evops(1,1)%Ufe(k,l, r,s, INT((i-1)/gt(1)) + 1) = conjg(rho(k,l,i))

					if((submethod == 's' .or. submethod == 'S') .and. (k /= r .or. l /= s)) then
						evops(1,1)%Ufe(k,l, r,s, INT((i-1)/gt(1)) + 1) = 0.0_dp
					end if
				endif
			end do

			write(cbuff,'(a,a1,a,i1,i1,i1,i1,a)') 'U_',type,'(',k,l,r,s,') filled'
					call print_log_message(trim(cbuff),5)


			end do
			end do

		end do
		end do

		call deinit_qme(submethod)

	end subroutine fill_evolution_superoperator_nakajima_zwanzig

	!*************************************************************
	!  Writing out Redfield tensor
	!*************************************************************

	subroutine write_redfield_tensor2(type)
		character, intent(in) :: type
		integer	(i4b)		:: i,j, t_index
		integer(i4b)		:: Uelement, Uelement2,Utnemele,Utnemele2, Ublock
		character(len=4)	:: no1,no2,no3,no4
		character(len=100)	:: name
		character(len=50)	:: prefix
		complex(dpc), dimension(:,:,:,:,:), pointer		:: actual_U
		complex(dpc), dimension(:,:,:,:,:), allocatable	:: Rdf
		complex(dpc), dimension(N1_from_type(type)*N2_from_type(type),N1_from_type(type)*N2_from_type(type))	:: Rdf_2ind

		Ublock = 1

		! We set indices range according to block we evaluate. Because rho0 is
		! whole density matrix, while evolution operators are only from particular
		! block, offset is set between these indices.
		if (type == '2') then
			actual_U => evops(Ublock,Ublock)%Ufe
			prefix = 'Redf2_fe'
		else if (type == 'E') then
			actual_U => evops(Ublock,Ublock)%Uee
			prefix = 'Redf2_ee'
		else if (type == 'O') then
			actual_U => evops(Ublock,Ublock)%Ueg
			prefix = 'Redf2_eg'
		end if

		ALLOCATE(Rdf,(size(actual_U,1),size(actual_U,2),size(actual_U,3),size(actual_U,4),size(actual_U,5) ) )

		do t_index=1, size(actual_U,5)
			call fill_superoperator_R(Rdf_2ind, t_index*gt(1), type)
			call superops_2indexed_to_4indexed(Rdf_2ind, Rdf(:,:,:,:,t_index), type)
		end do

		do Uelement=1,N1_from_type(type)
		do Uelement2=1,N2_from_type(type)
		do Utnemele=1,N1_from_type(type)
		do Utnemele2=1,N2_from_type(type)


		if(Uelement < 10) then
			write(no1,'(i1)')	Uelement
		else if (Uelement < 100) then
			write(no1,'(i2)')	Uelement
		else
			write(no1,'(i3)')	Uelement
		endif
		if(Uelement2 < 10) then
			write(no2,'(i1)')	Uelement2
		else if (Uelement2 < 100) then
			write(no2,'(i2)')	Uelement2
		else
			write(no2,'(i3)')	Uelement2
		endif
		if(Utnemele < 10) then
			write(no3,'(i1)')	Utnemele
		else if (Uelement2 < 100) then
			write(no3,'(i2)')	Utnemele
		else
			write(no3,'(i3)')	Utnemele
		endif
		if(Utnemele2 < 10) then
			write(no4,'(i1)')	Utnemele2
		else if (Uelement2 < 100) then
			write(no4,'(i2)')	Utnemele2
		else
			write(no4,'(i3)')	Utnemele2
		endif

		name = trim(prefix) // trim(no1) // '-'//trim(no2)//'--'// trim(no3) // '-'//trim(no4)//'.dat'

		open(UNIT=22, FILE = trim(file_join(out_dir,trim(name))))

		i = 1
		do while (i <= size(actual_U,5))
			write(22,*) gt(1)*dt*(i-1),' ',real(Rdf(Uelement,Uelement2,Utnemele,Utnemele2,i)),' ',aimag(Rdf(Uelement,Uelement2,Utnemele,Utnemele2,i))
			i = i + 1
		end do

		close(UNIT=22)

		end do
		end do
		end do
		end do

		DEALLOCATE(Rdf)
	end subroutine write_redfield_tensor2

end module qme_nakajima_zwanzig

