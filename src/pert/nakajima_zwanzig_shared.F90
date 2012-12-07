#include "util_allocation.h"

#define SUPERINDEX_FROM_K_L(k, l, lmax) ((l) + (lmax)*((k)-1))
#define K_FROM_SUPERINDEX(superindex, lmax) (((superindex) - 1) / (lmax) + 1)
#define L_FROM_SUPERINDEX(superindex, lmax) (mod(((superindex)-1) , (lmax)) + 1)

module nakajima_zwanzig_shared

	use std_types
	use std_io
	use numer_matrix
	use sci_misc

	use util_allocation
	use resources
	use numer_ode

	implicit none

	public::N1_from_type
	public::N2_from_type
	public::operator_to_exc
	public::operator_from_exc
	public::superops_to_exc ! transforms the superoperator into exc picture directly

	private::derivative_rates
	private::derivative_rates_inner
	public::prepare_very_specific_site_ops

	private::LQHO
	private::prepare_rho_exp
	private::clean_rho_exp
	private::rho_exp

	interface superops_to_exc
		module procedure superops_to_exc_2indexed
		module procedure superops_to_exc_4indexed
	end interface

	interface superops_from_exc
		module procedure superops_from_exc_2indexed
		module procedure superops_from_exc_4indexed
	end interface

	real(dp), dimension(:), pointer, private :: eng,en2

	complex(dpc), dimension(:,:), private, allocatable 		:: rhoexpXX, rhoexpxxx
	complex(dpc), dimension(:), private, allocatable		:: rhoexpYY
	real(dp), private :: k_et2, k_et1, k_ic21, k_sink, k_rate

	contains

    pure function N1_from_type(type) result(NN1)
    ! gives left operator size according to type
 		character, intent(in) 			:: type
 		integer(i4b)						:: NN1

 		if(type == 'g') then
 			NN1 = 1
 		else if(type == 'O') then
 			NN1 = iblocks(1,1)%eblock%N1
 		else if(type == 'E') then
 			NN1 = iblocks(1,1)%eblock%N1
 		else if(type == '2') then
 			NN1 = iblocks(1,1)%eblock%N1*(iblocks(1,1)%eblock%N1-1)/2
 		else if(type == 'F') then
 			NN1 = iblocks(1,1)%eblock%N1*(iblocks(1,1)%eblock%N1-1)/2
 		else
 			NN1 = -1
 		end if

 	end function

    pure function N2_from_type(type) result(NN2)
    ! gives right operator size according to type
 		character, intent(in) :: type
 		integer		:: NN2

 		if(type == 'g') then
 			NN2 = 1
 		else if(type == 'O') then
 			NN2 = 1
 		else if(type == 'E') then
 			NN2 = iblocks(1,1)%eblock%N1
 		else if(type == '2') then
 			NN2 = iblocks(1,1)%eblock%N1
 		else if(type == 'F') then
 			NN2 = iblocks(1,1)%eblock%N1*(iblocks(1,1)%eblock%N1-1)/2
 		else
 			NN2 = -1
 		end if

 	end function

 	subroutine operator_to_exc(fromto, type)
		complex(dpc), dimension(:,:), intent(inout) 				:: fromto
		character													:: type

		integer														:: NN1, NN2

		NN1 = N1_from_type(type)
		NN2 = N2_from_type(type)

		if(.not.(type == 'O' .or. type == 'E' .or. type == '2' .or. type == 'F')) then
 			call print_error_message(-1,  'wrong type in operator_to_exc() : ')
 			stop
 		else if(type == 'O') then
 			fromto = matmul(transpose(iblocks(1,1)%eblock%S1),fromto)
 		else if(type == 'E') then
 			fromto = matmul(matmul(iblocks(1,1)%eblock%S1,fromto),iblocks(1,1)%eblock%SS)
 		else if(type == '2') then
 			fromto = matmul(matmul(iblocks(1,1)%eblock%S1_2,fromto),iblocks(1,1)%eblock%SS)
 		else if(type == 'F') then
 			fromto = matmul(matmul(iblocks(1,1)%eblock%S1_2,fromto),iblocks(1,1)%eblock%SS_2)
 		end if

 	end subroutine operator_to_exc

 	subroutine operator_from_exc(fromto, type)
		complex(dpc), dimension(:,:), intent(inout) 				:: fromto
		character													:: type

		integer														:: NN1, NN2

		NN1 = N1_from_type(type)
		NN2 = N2_from_type(type)

		if(.not.(type == 'O' .or. type == 'E' .or. type == '2' .or. type == 'F')) then
 			call print_error_message(-1,  'wrong type in operator_to_exc() : ')
 			stop
 		else if(type == 'O') then
 			fromto = matmul(transpose(iblocks(1,1)%eblock%SS),fromto)
 		else if(type == 'E') then
 			fromto = matmul(matmul(iblocks(1,1)%eblock%SS,fromto),iblocks(1,1)%eblock%S1)
 		else if(type == '2') then
 			fromto = matmul(matmul(iblocks(1,1)%eblock%SS_2,fromto),iblocks(1,1)%eblock%S1)
 		else if(type == 'F') then
 			fromto = matmul(matmul(iblocks(1,1)%eblock%SS_2,fromto),iblocks(1,1)%eblock%S1_2)
 		end if

 	end subroutine operator_from_exc

	subroutine superops_to_exc_2indexed(fromto, type)
		complex(dpc), dimension(:,:), intent(inout)		:: fromto
		character											:: type

		complex(dpc), dimension(size(fromto,1),size(fromto,1)) :: Xi, Xi1
		integer											 	:: a,b,c,d
		integer												:: NN1, NN2
		complex(dpc), dimension(N1_from_type(type),N2_from_type(type)) :: rho

		NN1 = N1_from_type(type)
		NN2 = N2_from_type(type)
		Xi  = 0.0_dp
		Xi1 = 0.0_dp

		do d=1, NN2
		do c=1, NN1
				rho = 0.0_dp
				rho(c,d) = 1.0_dp

				call operator_to_exc(rho,type)

				do b=1, NN2
				do a=1, NN1
					Xi(SUPERINDEX_FROM_K_L(a,b,NN2),SUPERINDEX_FROM_K_L(c,d,NN2)) = rho(a,b)
				end do
				end do
		end do
		end do

		call inv(Xi,Xi1)

		fromto = matmul(matmul(Xi1,fromto),Xi)

	end subroutine superops_to_exc_2indexed

	subroutine superops_to_exc_4indexed(fromto, type)
		complex(dpc), dimension(:,:,:,:), intent(inout)	:: fromto
		character, intent(in)								:: type

		complex(dpc), dimension(N1_from_type(type)*N2_from_type(type), N1_from_type(type)*N2_from_type(type)) :: XX
		integer(i4b) :: i,j,k,l,  first_N, second_N

		XX = 0.0_dp

		first_N  = N1_from_type(type)
		second_N = N2_from_type(type)

		if(.not.(size(fromto,1) == first_N .or. size(fromto,2) == second_N .or. size(fromto,3) == first_N .or. size(fromto,4) == second_N)) then
 			call print_error_message(-1,  'wrong dimension in superops_to_exc_4indexed() : ')
 			stop
 		end if

		do k=1, first_N
		do l=1, second_N

		do i=1, first_N
		do j=1, second_N

		XX(SUPERINDEX_FROM_K_L(i,j,second_N),SUPERINDEX_FROM_K_L(k,l,second_N)) = &
			fromto(i,j,k,l)

		end do
		end do

		end do
		end do

		call superops_to_exc(XX, type)

		do k=1, first_N
		do l=1, second_N

		do i=1, first_N
		do j=1, second_N

		fromto(i,j,k,l) = XX(SUPERINDEX_FROM_K_L(i,j,second_N),SUPERINDEX_FROM_K_L(k,l,second_N))

		end do
		end do

		end do
		end do

	end subroutine superops_to_exc_4indexed

	subroutine superops_from_exc_2indexed(fromto, type)
		complex(dpc), dimension(:,:), intent(inout)		:: fromto
		character											:: type

		complex(dpc), dimension(size(fromto,1),size(fromto,1)) :: Xi, Xi1
		integer											 	:: a,b,c,d
		integer												:: NN1, NN2
		complex(dpc), dimension(N1_from_type(type),N2_from_type(type)) :: rho

		NN1 = N1_from_type(type)
		NN2 = N2_from_type(type)
		Xi  = 0.0_dp
		Xi1 = 0.0_dp

		do d=1, NN2
		do c=1, NN1
				rho = 0.0_dp
				rho(c,d) = 1.0_dp

				call operator_to_exc(rho,type)

				do b=1, NN2
				do a=1, NN1
					Xi(SUPERINDEX_FROM_K_L(a,b,NN2),SUPERINDEX_FROM_K_L(c,d,NN2)) = rho(a,b)
				end do
				end do
		end do
		end do

		call inv(Xi,Xi1)

		fromto = matmul(matmul(Xi,fromto),Xi1)

	end subroutine superops_from_exc_2indexed

	subroutine superops_from_exc_4indexed(fromto, type)
		complex(dpc), dimension(:,:,:,:), intent(inout)	:: fromto
		character, intent(in)								:: type

		complex(dpc), dimension(N1_from_type(type)*N2_from_type(type), N1_from_type(type)*N2_from_type(type)) :: XX
		integer(i4b) :: i,j,k,l,  first_N, second_N

		XX = 0.0_dp

		first_N  = N1_from_type(type)
		second_N = N2_from_type(type)

		if(.not.(size(fromto,1) == first_N .or. size(fromto,2) == second_N .or. size(fromto,3) == first_N .or. size(fromto,4) == second_N)) then
 			call print_error_message(-1,  'wrong dimension in superops_to_exc_4indexed() : ')
 			stop
 		end if

		do k=1, first_N
		do l=1, second_N

		do i=1, first_N
		do j=1, second_N

		XX(SUPERINDEX_FROM_K_L(i,j,second_N),SUPERINDEX_FROM_K_L(k,l,second_N)) = &
			fromto(i,j,k,l)

		end do
		end do

		end do
		end do

		call superops_from_exc(XX, type)

		do k=1, first_N
		do l=1, second_N

		do i=1, first_N
		do j=1, second_N

		fromto(i,j,k,l) = XX(SUPERINDEX_FROM_K_L(i,j,second_N),SUPERINDEX_FROM_K_L(k,l,second_N))

		end do
		end do

		end do
		end do

	end subroutine superops_from_exc_4indexed

	subroutine superops_4indexed_to_2indexed(from, to, type)
		complex(dpc), dimension(:,:,:,:), intent(in)		:: from
		character, intent(in)								:: type

		complex(dpc), dimension(:,:), intent(out) :: to
		integer(i4b) :: i,j,k,l,  first_N, second_N

		if(.not.(size(to,1) == size(to,2) .and. (N1_from_type(type)*N2_from_type(type) == size(to,1)) .and. &
			(size(from,1) == size(from,3)) .and. (size(from,2) == size(from,4)) .and. (size(from,1) == N1_from_type(type)) &
			.and. (size(from,2) == N2_from_type(type)) )) then

			call print_error_message(-1, "superops_4indexed_to_2indexed - size error")
			stop
		end if

		to = 0.0_dp

		first_N  = N1_from_type(type)
		second_N = N2_from_type(type)

		do k=1, first_N
		do l=1, second_N

		do i=1, first_N
		do j=1, second_N

		to(SUPERINDEX_FROM_K_L(i,j,second_N),SUPERINDEX_FROM_K_L(k,l,second_N)) = &
			from(i,j,k,l)

		end do
		end do

		end do
		end do

	end subroutine superops_4indexed_to_2indexed


	subroutine superops_2indexed_to_4indexed(from, to, type)
		complex(dpc), dimension(:,:), intent(in)		:: from
		character, intent(in)							:: type

		complex(dpc), dimension(:,:,:,:), intent(out) :: to
		integer(i4b) :: i,j,k,l,  first_N, second_N

		if(.not.(size(from,1) == size(from,2) .and. (N1_from_type(type)*N2_from_type(type) == size(from,1)) .and. &
			(size(to,1) == size(to,3)) .and. (size(to,2) == size(to,4)) .and. (size(to,1) == N1_from_type(type)) &
			.and. (size(to,2) == N2_from_type(type)) )) then

			call print_error_message(-1, "superops_2indexed_to_4indexed - size error")
			stop
		end if

		to = 0.0_dp

		first_N  = N1_from_type(type)
		second_N = N2_from_type(type)


		do k=1, first_N
		do l=1, second_N

		do i=1, first_N
		do j=1, second_N

		to(i,j,k,l) = from(SUPERINDEX_FROM_K_L(i,j,second_N),SUPERINDEX_FROM_K_L(k,l,second_N))

		end do
		end do

		end do
		end do

	end subroutine superops_2indexed_to_4indexed

	subroutine redfield_from_evops(from, to, type, tstep)
		complex(dpc), dimension(:,:,:,:,:), intent(in)	:: from
		complex(dpc), dimension(:,:,:,:,:), intent(out)	:: to
		character, intent(in)								:: type
		real(dp), intent(in)								:: tstep

		complex(dpc), dimension(N1_from_type(type)*N2_from_type(type), N1_from_type(type)*N2_from_type(type)) :: XX, YY, ZZ
		complex(dpc), dimension(:,:,:), allocatable :: xxx,yyy
		real(dp), dimension(:), allocatable :: tmatrix
		integer(i4b) :: i,j,k,l, NNN,  first_N, second_N

		XX = 0.0_dp

		first_N  = N1_from_type(type)
		second_N = N2_from_type(type)

		if(.not.(size(from,1) == first_N .or. size(from,2) == second_N .or. size(from,3) == first_N .or. size(from,4) == second_N &
				.or. size(from,1) == size(to,1) .or. size(from,2) == size(to,2) .or. size(from,3) == size(to,3) .or. size(from,4) == size(to,4) .or. size(from,5) == size(to,5) ) ) then
 			call print_error_message(-1,  'wrong dimension in superops_to_exc_4indexed() : ')
 			stop
 		end if

		ALLOCATE(xxx,(N1_from_type(type)*N2_from_type(type), N1_from_type(type)*N2_from_type(type),size(to,5)))
		ALLOCATE(yyy,(N1_from_type(type)*N2_from_type(type), N1_from_type(type)*N2_from_type(type),size(to,5)))
		ALLOCATE(tmatrix,(size(to,5)))

		do NNN = 1, size(to,5) ! do over time

		do k=1, first_N
		do l=1, second_N

		do i=1, first_N
		do j=1, second_N

		XX(SUPERINDEX_FROM_K_L(i,j,second_N),SUPERINDEX_FROM_K_L(k,l,second_N)) = &
			from(i,j,k,l,NNN)

		xxx(:,:,NNN) = XX(:,:)

		end do
		end do

		end do
		end do

		tmatrix(NNN) = (tstep)*(NNN-1)

		end do

		XX = (xxx(:,:, 2) - xxx(:,:,1))/(tstep)
		YY = 0.0_dp
		call spline(tmatrix,xxx,XX,YY,yyy)

		do NNN = 1, size(to,5) ! do over time

		! here is the procedure itself
		call splint_matrix(tmatrix,xxx,yyy,tmatrix(NNN),XX)
		call splint_matrix(tmatrix,xxx,yyy,tmatrix(NNN)+tstep/10.0_dp,YY)

		ZZ = -(XX - YY)/(tstep/10.0_dp)
		ZZ = (xxx(:,:, min(NNN+1,size(to,5)) ) - xxx(:,:,NNN))/(tstep)
		YY = transpose(conjg(XX))
		XX = matmul(ZZ,YY)

		do k=1, first_N
		do l=1, second_N

		do i=1, first_N
		do j=1, second_N

		to(i,j,k,l,NNN) = XX(SUPERINDEX_FROM_K_L(i,j,second_N),SUPERINDEX_FROM_K_L(k,l,second_N))

		end do
		end do

		end do
		end do

		end do ! do over time

		DEALLOCATE(xxx)
		DEALLOCATE(yyy)
		DEALLOCATE(tmatrix)

	end subroutine redfield_from_evops

	!*************************************************************
	!  Writing out evolution operators
	!*************************************************************

	subroutine write_evolution_operators(type)
		character, intent(in) :: type
		integer	(i4b)		:: i,j
		integer(i4b)		:: Uelement, Uelement2,Utnemele,Utnemele2, Ublock
		character(len=4)	:: no1,no2,no3,no4
		character(len=100)	:: name
		character(len=50)	:: prefix
		complex(dpc), dimension(:,:,:,:,:), pointer		:: actual_U

		Ublock = 1

		! We set indices range according to block we evaluate. Because rho0 is
		! whole density matrix, while evolution operators are only from particular
		! block, offset is set between these indices.
		if (type == '2') then
			actual_U => evops(Ublock,Ublock)%Ufe
			prefix = 'Evops_fe'
		else if (type == 'E') then
			actual_U => evops(Ublock,Ublock)%Uee
			prefix = 'Evops_ee'
		else if (type == 'O') then
			actual_U => evops(Ublock,Ublock)%Ueg
			prefix = 'Evops_eg'
		end if

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
			write(22,*) gt(1)*dt*(i-1),' ',real(actual_U(Uelement,Uelement2,Utnemele,Utnemele2,i)),' ',aimag(actual_U(Uelement,Uelement2,Utnemele,Utnemele2,i))
			i = i + 1
		end do

		close(UNIT=22)

		end do
		end do
		end do
		end do
	end subroutine write_evolution_operators


	!*************************************************************
	!  Writing out Redfield tensor
	!*************************************************************

	subroutine write_redfield_tensor(type)
		character, intent(in) :: type
		integer	(i4b)		:: i,j
		integer(i4b)		:: Uelement, Uelement2,Utnemele,Utnemele2, Ublock
		character(len=4)	:: no1,no2,no3,no4
		character(len=100)	:: name
		character(len=50)	:: prefix
		complex(dpc), dimension(:,:,:,:,:), pointer		:: actual_U
		complex(dpc), dimension(:,:,:,:,:), allocatable	:: Rdf

		Ublock = 1

		! We set indices range according to block we evaluate. Because rho0 is
		! whole density matrix, while evolution operators are only from particular
		! block, offset is set between these indices.
		if (type == '2') then
			actual_U => evops(Ublock,Ublock)%Ufe
			prefix = 'Redf_fe'
		else if (type == 'E') then
			actual_U => evops(Ublock,Ublock)%Uee
			prefix = 'Redf_ee'
		else if (type == 'O') then
			actual_U => evops(Ublock,Ublock)%Ueg
			prefix = 'Redf_eg'
		end if

		ALLOCATE(Rdf,(size(actual_U,1),size(actual_U,2),size(actual_U,3),size(actual_U,4),size(actual_U,5) ) )

		call redfield_from_evops(actual_U, Rdf, type, gt(1)*dt)

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
	end subroutine write_redfield_tensor

	! vibrational levels
	integer(i4b) function to_vibrational_multiindex(levels_list) result(multiindex)
		integer(i4b), dimension(:), intent(in) :: levels_list
		integer(i4b) :: i

		if(size(levels_list) /= size(current_s_block%QHO_lvls)) then
			call print_error_message(-1, "dimension error in to_vibrational_multiindex")
		end if
		if(maxval(levels_list - current_s_block%QHO_lvls) >= 0 .or. minval(levels_list) < 0) then
			call print_error_message(-1, "dimension error in to_vibrational_multiindex")
		end if

		multiindex = 0
		do i=1, size(levels_list)
			multiindex = multiindex + levels_list(i)
			multiindex = multiindex * current_s_block%QHO_lvls(i)
		end do
		multiindex = multiindex / current_s_block%QHO_lvls(i-1) + 1

	end function to_vibrational_multiindex

	function from_vibrational_multiindex(multiindex) result(levels_list)
		integer(i4b), dimension(size(current_s_block%QHO_lvls)) :: levels_list
		integer(i4b), intent(in) :: multiindex
		integer(i4b) :: i,n

		levels_list = 0
		n = multiindex - 1
		n = n * current_s_block%QHO_lvls(size(levels_list))
		do i=size(levels_list),1,-1
			n = (n - levels_list(i))/current_s_block%QHO_lvls(i)
			levels_list(i) = mod(n, current_s_block%QHO_lvls(i))
		end do

	end function from_vibrational_multiindex

	subroutine add_vibrational_DOF(X,Y,type)
		character, intent(in) :: type
		real(dp), dimension(:,:), intent(in) :: X

		real(dp), dimension(:,:), intent(out) :: Y
		integer	(i4b)		:: i,j,a,b,t,u,r,s, N0vib
		integer	(i4b), dimension(size(current_s_block%QHO_lvls))	:: vi1,vi2,vvi1,vvi2
		real(dp) :: mult, d1, d2, summ

		write(*,*) size(X,1), size(X,2), size(Y,1), size(Y,2)
		write(*,*) N1_from_type(type), N2_from_type(type), N1_from_type(type)*product(current_s_block%QHO_lvls), N2_from_type(type)*product(current_s_block%QHO_lvls)

		! PARAMETERS OF FC FACTORS CALCULATION NOT TESTED !

		N0vib = product(current_s_block%QHO_lvls)

		if(size(X,1) /= N1_from_type(type) .or.  size(X,2) /= N2_from_type(type) .or. &
			size(Y,1) /= N1_from_type(type)*N0vib .or.  &
			size(Y,2) /= N2_from_type(type)*N0vib) then
			call print_error_message(-1, "dimension error in add_vibrational_DOF")
		end if

		if(type == 'g') then

			Y = 0.0_dp
			do a=1, N1_from_type(type)
			do b=1, N2_from_type(type)
			do i=1, N0vib

			j = i

			Y(i+(a-1)*N0vib,j+(b-1)*N0vib) = X(a,b)

			end do
			end do
			end do

		else if(type == 'E' .or. type == 'O') then
			do a=1, N1_from_type(type)
			do b=1, N2_from_type(type)
			do i=1, N0vib
			do j=1, N0vib

				vi1 = from_vibrational_multiindex(i)
				vi2 = from_vibrational_multiindex(j)

				summ = 0.0_dp
				do r=1, N0vib
!				do s=1, N0vib
!
!				if(r /= s) then
!					cycle
!				end if
				s = r

				vvi1 = from_vibrational_multiindex(r)
				vvi2 = from_vibrational_multiindex(s)

				mult = 1.0_dp
				do t=1,size(vi1)

				d1 = sqrt(2*current_s_block%QHO_hrfact(t))
				d2 = sqrt(2*current_s_block%QHO_hrfact(t))

				if(vi2(t) /= vvi2(t) .and. type == 'O') then
					mult = 0.0_dp
				end if

				mult = mult * franc_condon_factor(vi1(t),current_s_block%QHO_freq(t),vvi1(t),current_s_block%QHO_freq(t),d1)

!				write(*,*) '*  ', vi1(t),current_s_block%QHO_freq(t),vvi1(t),current_s_block%QHO_freq(t),franc_condon_factor(vi1(t),current_s_block%QHO_freq(t),vvi1(t),current_s_block%QHO_freq(t),d1)
!				write(*,*) '*  ', vi2(t),current_s_block%QHO_freq(t),vvi2(t),current_s_block%QHO_freq(t),franc_condon_factor(vi2(t),current_s_block%QHO_freq(t),vvi2(t),current_s_block%QHO_freq(t),d2)
!				write(*,*)

				end do

!				write(*,*) mult
!				write(*,*) '---------------------------'
				summ = summ + mult

!				end do
				end do

!				write(*,*) 'summ=',summ
!				write(*,*) '---------------------------'
				Y(i+(a-1)*N0vib,j+(b-1)*N0vib) = X(a,b)*summ


			end do
			end do
			end do
			end do
		else if(type == 'F') then
			call print_error_message(-1,"I'm not prepared for this!")
		end if


	end subroutine add_vibrational_DOF

	real(dp) function LQHO(in,what,in2) result(nav)
		integer(i4b), intent(in) :: in
		character(len=*), intent(in):: what
		integer(i4b), intent(in), optional :: in2
		real(dp), parameter :: &
			LQHO_freq_LO = 1150.0_dp/Energy_internal_to_cm, &
			LQHO_freq_HI = 1520.0_dp/Energy_internal_to_cm
		real(dp) ::	LQHO_HR_LO = 0.6_dp, &
					LQHO_HR_HI = 0.6_dp
		real(dp) :: d_LO,d_HI

		! Fixed vibrational structure, parameters input here
		!    LO HI
		!  1 |0>|0>
		!  2 |1>|0>
		!  3 |0>|1>
		!  4 |2>|0>
		!  5 |1>|1>
		!  6 |0>|2>

		! possible transitions:
		!   1-1
		!   2-2
		!   3-3
		!   4-4
		!   5-5
		!   6-6

		!   1-2
		!   1-3
		!   2-4
		!   2-5
		!   3-6
		!   3-5

		if(what /= "EN" .and. what /= "DD" .and. what /= "D2") then
			call print_error_message(-1,"Wrong 'what' variable in LHQO(...)")
		end if

		if(what == 'DD' .or. what == 'EN') then
			LQHO_HR_LO = current_s_block%QHO_hrfact(1)
			LQHO_HR_HI = current_s_block%QHO_hrfact(4)
		elseif(what == 'D2') then
			LQHO_HR_LO = 0.2
			LQHO_HR_HI = 0.3
		end if

		nav = 0.0_dp
		d_LO = sqrt(2*LQHO_HR_LO/LQHO_freq_LO)
		d_HI = sqrt(2*LQHO_HR_HI/LQHO_freq_HI)

		if(what == 'EN') then

			if(in == 1) then
				nav = 0.0_dp
			elseif(in == 2) then
				nav = LQHO_freq_LO
			elseif(in == 3) then
				nav = LQHO_freq_HI
			elseif(in == 4) then
				nav = LQHO_freq_LO*2
			elseif(in == 5) then
				nav = LQHO_freq_LO + LQHO_freq_HI
			elseif(in == 6) then
				nav = LQHO_freq_HI*2
			end if

		elseif(what == 'DD' .or. what == 'D2') then
			if(.not.present(in2)) then
				call print_error_message(-1, "third argument must be given for 'DD' in  LQHO()")
			end if

			if(in == 1 .and. in2 == 1) then
				nav = franc_condon_factor(0,LQHO_freq_LO,0,LQHO_freq_LO,d_LO)*&
						franc_condon_factor(0,LQHO_freq_HI,0,LQHO_freq_HI,d_HI)
			elseif(in == 2 .and. in2 == 2) then
				nav = franc_condon_factor(1,LQHO_freq_LO,1,LQHO_freq_LO,d_LO)*&
						franc_condon_factor(0,LQHO_freq_HI,0,LQHO_freq_HI,d_HI)           *1
			elseif(in == 3 .and. in2 == 3) then
				nav = franc_condon_factor(0,LQHO_freq_LO,0,LQHO_freq_LO,d_LO)*&
						franc_condon_factor(1,LQHO_freq_HI,1,LQHO_freq_HI,d_HI)                              *1
			elseif(in == 4 .and. in2 == 4) then
				nav = franc_condon_factor(2,LQHO_freq_LO,2,LQHO_freq_LO,d_LO)*&
						franc_condon_factor(0,LQHO_freq_HI,0,LQHO_freq_HI,d_HI)           *1
			elseif(in == 5 .and. in2 == 5) then
				nav = franc_condon_factor(1,LQHO_freq_LO,1,LQHO_freq_LO,d_LO)*&
						franc_condon_factor(1,LQHO_freq_HI,1,LQHO_freq_HI,d_HI)           *1                 *1
			elseif(in == 6 .and. in2 == 6) then
				nav = franc_condon_factor(0,LQHO_freq_LO,0,LQHO_freq_LO,d_LO)*&
						franc_condon_factor(2,LQHO_freq_HI,2,LQHO_freq_HI,d_HI)                              *1


			elseif((in == 1 .and. in2 == 2) .or. (in == 2 .and. in2 == 1)) then
				nav = franc_condon_factor(  0,LQHO_freq_LO,1,LQHO_freq_LO,d_LO)*&
						franc_condon_factor(0,LQHO_freq_HI,0,LQHO_freq_HI,d_HI)           *1
			elseif((in == 1 .and. in2 == 3) .or. (in == 3 .and. in2 == 1)) then
				nav = franc_condon_factor(  0,LQHO_freq_LO,0,LQHO_freq_LO,d_LO)*&
						franc_condon_factor(0,LQHO_freq_HI,1,LQHO_freq_HI,d_HI)                              *1
			elseif((in == 1 .and. in2 == 4) .or. (in == 4 .and. in2 == 1)) then
				nav = franc_condon_factor(  0,LQHO_freq_LO,2,LQHO_freq_LO,d_LO)*&
						franc_condon_factor(0,LQHO_freq_HI,0,LQHO_freq_HI,d_HI)           *1
			elseif((in == 1 .and. in2 == 5) .or. (in == 5 .and. in2 == 1)) then
				nav = franc_condon_factor(  0,LQHO_freq_LO,1,LQHO_freq_LO,d_LO)*&
						franc_condon_factor(0,LQHO_freq_HI,1,LQHO_freq_HI,d_HI)           *1                 *1
			elseif((in == 1 .and. in2 == 6) .or. (in == 6 .and. in2 == 1)) then
				nav = franc_condon_factor(  0,LQHO_freq_LO,0,LQHO_freq_LO,d_LO)*&
						franc_condon_factor(0,LQHO_freq_HI,2,LQHO_freq_HI,d_HI)                              *1

			elseif((in == 3 .and. in2 == 2) .or. (in == 2 .and. in2 == 3)) then
				nav = franc_condon_factor(  1,LQHO_freq_LO,0,LQHO_freq_LO,d_LO)*&
						franc_condon_factor(0,LQHO_freq_HI,1,LQHO_freq_HI,d_HI)           *1                 *1
			elseif((in == 4 .and. in2 == 2) .or. (in == 2 .and. in2 == 4)) then
				nav = franc_condon_factor(  1,LQHO_freq_LO,2,LQHO_freq_LO,d_LO)*&
						franc_condon_factor(0,LQHO_freq_HI,0,LQHO_freq_HI,d_HI)           *1
			elseif((in == 5 .and. in2 == 2) .or. (in == 2 .and. in2 == 5)) then
				nav = franc_condon_factor(  1,LQHO_freq_LO,1,LQHO_freq_LO,d_LO)*&
						franc_condon_factor(0,LQHO_freq_HI,1,LQHO_freq_HI,d_HI)           *1                 *1
			elseif((in == 6 .and. in2 == 2) .or. (in == 2 .and. in2 == 6)) then
				nav = franc_condon_factor(  1,LQHO_freq_LO,0,LQHO_freq_LO,d_LO)*&
						franc_condon_factor(0,LQHO_freq_HI,2,LQHO_freq_HI,d_HI)           *1                 *1

			elseif((in == 4 .and. in2 == 3) .or. (in == 3 .and. in2 == 4)) then
				nav = franc_condon_factor(  0,LQHO_freq_LO,2,LQHO_freq_LO,d_LO)*&
						franc_condon_factor(1,LQHO_freq_HI,0,LQHO_freq_HI,d_HI)           *1                 *1
			elseif((in == 5 .and. in2 == 3) .or. (in == 3 .and. in2 == 5)) then
				nav = franc_condon_factor(  0,LQHO_freq_LO,1,LQHO_freq_LO,d_LO)*&
						franc_condon_factor(1,LQHO_freq_HI,1,LQHO_freq_HI,d_HI)           *1                 *1
			elseif((in == 6 .and. in2 == 3) .or. (in == 3 .and. in2 == 6)) then
				nav = franc_condon_factor(  0,LQHO_freq_LO,0,LQHO_freq_LO,d_LO)*&
						franc_condon_factor(1,LQHO_freq_HI,2,LQHO_freq_HI,d_HI)                              *1

			elseif((in == 4 .and. in2 == 5) .or. (in == 5 .and. in2 == 4)) then
				nav = franc_condon_factor(  2,LQHO_freq_LO,1,LQHO_freq_LO,d_LO)*&
						franc_condon_factor(0,LQHO_freq_HI,1,LQHO_freq_HI,d_HI)           *1                 *1
			elseif((in == 4 .and. in2 == 6) .or. (in == 6 .and. in2 == 4)) then
				nav = franc_condon_factor(  2,LQHO_freq_LO,0,LQHO_freq_LO,d_LO)*&
						franc_condon_factor(0,LQHO_freq_HI,2,LQHO_freq_HI,d_HI)           *1                 *1

			elseif((in == 6 .and. in2 == 5) .or. (in == 5 .and. in2 == 6)) then
				nav = franc_condon_factor(  1,LQHO_freq_LO,0,LQHO_freq_LO,d_LO)*&
						franc_condon_factor(1,LQHO_freq_HI,2,LQHO_freq_HI,d_HI)           *1                 *1

			end if
		end if
	end function LQHO

	subroutine prepare_very_specific_site_ops()
		real(dp), dimension(:,:), allocatable :: dd,dd2,J,dx,dx2,dy,dy2,dz,dz2
		real(dp), dimension(:), pointer :: eng_,en_,en2_
		integer	(i4b)		:: i,a,b,t,u,r,s,c,d,aa,bb,t_index, Nvvg, Nvve, Nvvf, Ncar
		real(dp)			:: d1, S1_en, tt
		complex(dpc), allocatable, dimension(:,:)		:: rho0, drho
		complex(dpc), allocatable, dimension(:,:,:)		:: rho
		character(len=256) :: cbuff
		complex(dpc), pointer, dimension(:) :: gg_tr1, gg_tr2, gg_tr3, gg_tr4, gg_tr5, gg_tr6
		logical :: evops_runge_kutta

		if(size(current_s_block%QHO_lvls) /= 6) then
			call print_error_message(-1, "prepare_very_specific_site_ops, wrong number of transitions (!= 5)")
		end if

		if(current_s_block%QHO_lvls(1) /= current_s_block%QHO_lvls(2) .or. current_s_block%QHO_lvls(2) /= current_s_block%QHO_lvls(3)) then
			call print_error_message(-1, "transition 1, 2 and 3 have to have the same number of levels")
		end if

		!
		! Ne block:
		!    S2
		!    S1
		!    Qx
		!    Qy
		!
		! Nf block
		!    Sn
		!    Qx-S1
		!    Qx-S2
		!    Qy-S2
		!    Sm
		!
		! Transitions
		!    1  g-S2
		!    2  g-S1
		!    3  S1-Sn
		!    4  g-Qx
		!    5  S2-Sm
		!	 6  g-Qy
		!
		! Fixed vibrational structure, parameters input here
		!    LO HI
		!    |0>|0>
		!    |1>|0>
		!    |0>|1>
		!    |2>|0>
		!    |1>|1>
		!    |0>|2>

		!Ncar = current_s_block%QHO_lvls(1)
		Ncar = 6 ! LHQO

		! only vibrational ground states
		Nvvg = Ncar

		! excited vstates of first transition (S2) + ground vstates of second transition (S1) + 2 Chl
		Nvve = Ncar*4

		! excited vstates of second transition (Sn) + ???
		Nvvf = Ncar*5

		special_carotenoid_hack = .true.
		DEALLOCATE(current_e_block%ione1)
		DEALLOCATE(current_e_block%ione2)
		ALLOCATE(current_e_block%ione1,(1)) ! is not used, will trigger alarm if is
		ALLOCATE(current_e_block%ione2,(1))

		ALLOCATE(dd,(Nvve,Nvvg))
		ALLOCATE(dd2,(Nvvf,Nvve))
		ALLOCATE(dx,(Nvve,Nvvg))
		ALLOCATE(dx2,(Nvvf,Nvve))
		ALLOCATE(dy,(Nvve,Nvvg))
		ALLOCATE(dy2,(Nvvf,Nvve))
		ALLOCATE(dz,(Nvve,Nvvg))
		ALLOCATE(dz2,(Nvvf,Nvve))
		ALLOCATE(J,(Nvve,Nvve))
		ALLOCATE(eng_,(Nvvg))
		ALLOCATE(en_,(Nvve))
		ALLOCATE(en2_,(Nvvf))

		en => en_
		eng => eng_
		en2 => en2_

		J = 0.0_dp
		S1_en = 0.0_dp

		dd = 0.0_dp

		do a=1, Nvve
			if(a <= Ncar) then
				en(a) = LQHO(a,"EN") + current_s_block%en(1) ! S2
			elseif(a <= 2*Ncar) then
				en(a) = LQHO(a-Ncar,"EN") + current_s_block%en(2) ! S1
			elseif(a <= 3*Ncar) then
				en(a) = LQHO(a-2*Ncar,"EN") + current_s_block%en(4) ! Qx Chl
			elseif(a <= 4*Ncar) then
				en(a) = LQHO(a-3*Ncar,"EN") + current_s_block%en(6) ! Qy Chl
			end if

		do b=1, Nvvg

			!eng(b) = (b-1)*current_s_block%QHO_freq(1)
			eng(b) = LQHO(b,"EN")

			if(a <= Ncar) then ! g - S2
				d1 = sqrt(2*current_s_block%QHO_hrfact(1)/current_s_block%QHO_freq(1))
				!dd(a,b) = current_s_block%dd(1,1)*franc_condon_factor(a-1,current_s_block%QHO_freq(1),b-1,current_s_block%QHO_freq(1),d1)
				!dx(a,b) = current_s_block%dx(1,1)*franc_condon_factor(a-1,current_s_block%QHO_freq(1),b-1,current_s_block%QHO_freq(1),d1)
				!dy(a,b) = current_s_block%dy(1,1)*franc_condon_factor(a-1,current_s_block%QHO_freq(1),b-1,current_s_block%QHO_freq(1),d1)
				!dz(a,b) = current_s_block%dz(1,1)*franc_condon_factor(a-1,current_s_block%QHO_freq(1),b-1,current_s_block%QHO_freq(1),d1)
				dd(a,b) = current_s_block%dd(1,1)*LQHO(a,"DD",b)
				dx(a,b) = current_s_block%dx(1,1)*LQHO(a,"DD",b)
				dy(a,b) = current_s_block%dy(1,1)*LQHO(a,"DD",b)
				dz(a,b) = current_s_block%dz(1,1)*LQHO(a,"DD",b)
			elseif(a <= 2*Ncar) then ! g - S1
				dd(a,b) = 0.0_dp
				dx(a,b) = 0.0_dp
				dy(a,b) = 0.0_dp
				dz(a,b) = 0.0_dp
			elseif(a <= 3*Ncar) then ! Qx
				aa = a - 2*Ncar
				d1 = 0 ! ground state to ground state

				if(aa == b) then	! HR = 0
					dd(a,b) = current_s_block%dd(4,1)*1
					dx(a,b) = current_s_block%dx(4,1)*1
					dy(a,b) = current_s_block%dy(4,1)*1
					dz(a,b) = current_s_block%dz(4,1)*1
				else
					dd(a,b) = current_s_block%dd(4,1)*0
					dx(a,b) = current_s_block%dx(4,1)*0
					dy(a,b) = current_s_block%dy(4,1)*0
					dz(a,b) = current_s_block%dz(4,1)*0
				end if
			elseif(a <= 4*Ncar) then ! g - Qy
				dd(a,b) = 0.0_dp
				dx(a,b) = 0.0_dp
				dy(a,b) = 0.0_dp
				dz(a,b) = 0.0_dp
			end if

		end do
		end do

		do a=1, Nvvf
			if(a <= Ncar) then
				aa = a
				en2(a) = LQHO(aa,"EN") + current_s_block%en(2) + current_s_block%en(3)
			elseif(a <= 2*Ncar) then
				aa = a - Ncar
				en2(a) = LQHO(aa,"EN") + current_s_block%en(2) + current_s_block%en(4)
			elseif(a <= 3*Ncar) then
				aa = a - 2*Ncar
				en2(a) = LQHO(aa,"EN") + current_s_block%en(1) + current_s_block%en(4)
			elseif(a <= 4*Ncar) then
				aa = a - 3*Ncar
				en2(a) = LQHO(aa,"EN") + current_s_block%en(1) + current_s_block%en(6)
			elseif(a <= 5*Ncar) then
				aa = a - 4*Ncar
				en2(a) = LQHO(aa,"EN") + current_s_block%en(1) + current_s_block%en(5)
			end if

		do b=1, Nvve

			if(b > Ncar .and. b <= 2*Ncar .and. a <= Ncar) then ! S1 - Sn
				d1 = sqrt(2*current_s_block%QHO_hrfact(3)/current_s_block%QHO_freq(3))
				aa = a
				bb = b - Ncar
				dd2(a,b) = current_s_block%dd(3,1)*LQHO(aa,"D2",bb)
				dx2(a,b) = current_s_block%dx(3,1)*LQHO(aa,"D2",bb)
				dy2(a,b) = current_s_block%dy(3,1)*LQHO(aa,"D2",bb)
				dz2(a,b) = current_s_block%dz(3,1)*LQHO(aa,"D2",bb)

			elseif(a > 4*Ncar .and. a <= 5*Ncar .and. b <= Ncar) then ! S2 - Sm
				d1 = sqrt(2*current_s_block%QHO_hrfact(5)/current_s_block%QHO_freq(5))
				aa = a - 4*Ncar
				bb = b
				dd2(a,b) = current_s_block%dd(5,1)*LQHO(aa,"D2",bb)
				dx2(a,b) = current_s_block%dx(5,1)*LQHO(aa,"D2",bb)
				dy2(a,b) = current_s_block%dy(5,1)*LQHO(aa,"D2",bb)
				dz2(a,b) = current_s_block%dz(5,1)*LQHO(aa,"D2",bb)
			elseif(a > 2*Ncar .and. a <= 3*Ncar .and. b > 2*Ncar .and. b <= 3*Ncar) then ! QxS2-Qx
				dd2(a,b) = dd(a-2*Ncar,b-2*Ncar)
				dx2(a,b) = dx(a-2*Ncar,b-2*Ncar)
				dy2(a,b) = dy(a-2*Ncar,b-2*Ncar)
				dz2(a,b) = dz(a-2*Ncar,b-2*Ncar)
			elseif(a > 2*Ncar .and. a <= 3*Ncar .and. b <= Ncar) then ! QxS2-S2
				dd2(a,b) = dd(a,b)
				dx2(a,b) = dx(a,b)
				dy2(a,b) = dy(a,b)
				dz2(a,b) = dz(a,b)
			elseif(a > Ncar .and. a <= 2*Ncar .and. b > 2*Ncar .and. b <= 3*Ncar) then ! QxS1-Qx
				dd2(a,b) = dd(a,b-2*Ncar)
				dx2(a,b) = dx(a,b-2*Ncar)
				dy2(a,b) = dy(a,b-2*Ncar)
				dz2(a,b) = dz(a,b-2*Ncar)
			elseif(a > Ncar .and. a <= 2*Ncar .and. b > Ncar .and. b <= 2*Ncar) then ! QxS1-S1
				dd2(a,b) = dd(a+Ncar,b-Ncar)
				dx2(a,b) = dx(a+Ncar,b-Ncar)
				dy2(a,b) = dy(a+Ncar,b-Ncar)
				dz2(a,b) = dz(a+Ncar,b-Ncar)
			elseif(a > 3*Ncar .and. a <= 4*Ncar .and. b > 3*Ncar .and. b <= 4*Ncar) then ! QyS2-Qy
				dd2(a,b) = dd(a-3*Ncar,b-3*Ncar)
				dx2(a,b) = dx(a-3*Ncar,b-3*Ncar)
				dy2(a,b) = dy(a-3*Ncar,b-3*Ncar)
				dz2(a,b) = dz(a-3*Ncar,b-3*Ncar)

				! fill in the rest !!
			else
				dd2(a,b) = 0.0_dp
				dx2(a,b) = 0.0_dp
				dy2(a,b) = 0.0_dp
				dz2(a,b) = 0.0_dp
			end if

		end do
		end do

		write(*,*) 'HRf',current_s_block%QHO_hrfact
		write(*,*)
		write(*,*) 'freq',current_s_block%QHO_freq
		write(*,*)
		write(*,*) 'dd',dd
		write(*,*)
		write(*,*) 'dd2',dd2
		write(*,*)
		write(*,*) 'eng',eng
		write(*,*)
		write(*,*) 'en', en
		write(*,*)
		write(*,*) 'en2', en2
		write(*,*)

		write(*,*) 'allocated? Ufe Ufes', associated(evops(1,1)%Ufe),associated(evops(1,1)%UfeS)


		! now, we multiply the dipole moments by the pulse spectrum
		do a=1, Nvve
		do b=1, Nvvg
			dd(a,b) = dd(a,b) * pulse_spectrum(en(a) - eng(b))
		end do
		end do

		do a=1, Nvvf
		do b=1, Nvve
			dd2(a,b) = dd2(a,b) * pulse_spectrum(en2(a) - en(b))
		end do
		end do

		write(*,*)
		write(*,*) 'dd',dd
		write(*,*)
		write(*,*) 'dd2',dd2



! ----------------------------- EVOPS -------------------------



		DEALLOCATE(evops(1,1)%Ueg)
		DEALLOCATE(evops(1,1)%Uee)
		DEALLOCATE(evops(1,1)%Ufe)
		DEALLOCATE(evops(1,1)%UfeS)
		DEALLOCATE(evops(1,1)%Ugg)

		ALLOCATE(evops(1,1)%Ueg,(Nvve, Nvvg,Nvve, Nvvg,Nt(1)))
		ALLOCATE(evops(1,1)%Uee,(Nvve,Nvve,Nvve,Nvve,Nt(2)))
        !ALLOCATE(evops(1,1)%Ufe,(Nvvf,Nvve,Nvvf,Nvve,Nt(1)))
 		ALLOCATE(evops(1,1)%UfeS,(Nvvf,Nvve,Nt(1)))
		ALLOCATE(evops(1,1)%Ugg,(Nvvg,Nvvg,Nvvg,Nvvg,Nt(2)))

	    evops(1,1)%Ueg = 0.0_dp
	    evops(1,1)%Uee = 0.0_dp
	    evops(1,1)%Ugg = 0.0_dp
	    evops(1,1)%UfeS = 0.0_dp
	    !evops(1,1)%Ufe = 0.0_dp

		gg_tr1 => igofts(iblocks(1,1)%sblock%gindex(1))%goft%gt
		gg_tr2 => igofts(iblocks(1,1)%sblock%gindex(2))%goft%gt
		gg_tr3 => igofts(iblocks(1,1)%sblock%gindex(3))%goft%gt
		gg_tr4 => igofts(iblocks(1,1)%sblock%gindex(4))%goft%gt
		gg_tr5 => igofts(iblocks(1,1)%sblock%gindex(5))%goft%gt
		gg_tr6 => igofts(iblocks(1,1)%sblock%gindex(6))%goft%gt


	    !coherences
	    do a=1, size(evops(1,1)%Ugg,1)
	    do b=1, size(evops(1,1)%Ugg,2)
	    !do c=1, size(evops(1,1)%Ugg,3)
	    !do d=1, size(evops(1,1)%Ugg,4)
	    if(a == b) then
	    	cycle
	    end if
	    c = a
	    d = b
	    do t=1, size(evops(1,1)%Ugg,5)

	    	evops(1,1)%Ugg(a,b,c,d,t) = exp(-(0,1)*(t-1)*dt*gt(1) &
	    		* (eng(a) - eng(b)) )

	    end do
	    !end do
	    !end do
	    end do
	    end do

		!populations
	    do a=1, size(evops(1,1)%Ugg,1)
	    !do b=1, size(evops(1,1)%Ugg,2)
	    do c=1, size(evops(1,1)%Ugg,3)
	    !do d=1, size(evops(1,1)%Ugg,4)
	    b = a
	    d = c
	    do t=1, size(evops(1,1)%Ugg,5)

	    if(a == c) then
	    	evops(1,1)%Ugg(a,b,c,d,t) = 1.0_dp
	    end if

	    end do
	    !end do
	    !end do
	    end do
	    end do

	    do a=1, size(evops(1,1)%Ueg,1)
	    do b=1, size(evops(1,1)%Ueg,2)
	    !do c=1, size(evops(1,1)%Ueg,3)
	    !do d=1, size(evops(1,1)%Ueg,4)
	    c = a
	    d = b
	    do t=1, size(evops(1,1)%Ueg,5)

	    	evops(1,1)%Ueg(a,b,c,d,t) = exp(-(0,1)*(t-1)*dt*gt(1) &
	    		* (en(a) - eng(b) - rwa) )

	    	if(a <= Ncar) then
	    		evops(1,1)%Ueg(a,b,c,d,t) = evops(1,1)%Ueg(a,b,c,d,t) * exp(-gg_tr1((t-1)*gt(1)+1))
	    	else if(a >= 2*Ncar+1 .and. a <= 3*Ncar) then
	    		evops(1,1)%Ueg(a,b,c,d,t) = evops(1,1)%Ueg(a,b,c,d,t) * exp(-gg_tr4((t-1)*gt(1)+1))
	    	end if

	    !end do
	    !end do
	    end do
	    end do
	    end do

		if(Nt(2) > 20) then
			evops_runge_kutta = .true.
		else
			evops_runge_kutta = .false.
		end if


		write(*,*) 'grid_Nt', grid_Nt
		write(*,*) 'evops_runge_kutta', evops_runge_kutta
		call flush()

		ALLOCATE(rho,(Nvve,Nvve,2))
		ALLOCATE(rho0,(Nvve,Nvve))
		ALLOCATE(drho,(Nvve,Nvve))

		call read_derivative_rates_from_file()

		if(.not. evops_runge_kutta) then
			call prepare_rho_exp(rho(:,:,1))
		end if

		evops(1,1)%Uee = 0.0_dp


		do a=1,Nvve
		do b=1,Nvve

		evops(1,1)%Uee(a,a,a,a,1) = 1.0_dp

		write(cbuff,'(a,i2,a,i2,a)') 'Evaluating Evops  ',a,' ',b,', :,  :'
		call print_log_message(trim(cbuff),5)


		rho = 0.0_dp
		rho(a,b,1) = 1.0_dp
		drho = 0.0_dp

		if (resources_output_contains("2d_ftpe")) then
		do t_index = 1, grid_Nt-1
			tt = (t_index-1)*dt

			if(evops_runge_kutta) then
				! Matrix Runge-Kutta -- slow but general and works well
				call derivative_rates(tt,rho(:,:,1),drho)
				call ode_rk4(rho(:,:,1),drho,tt,dt,rho(:,:,1+1),derivative_rates)
				rho(:,:,1) = rho(:,:,1+1)
			end if

			if(mod(t_index,gt(2)) == 0 .and. t_index/gt(2)+1 <= size(evops(1,1)%Uee,5)) then
				if(.not. evops_runge_kutta) then
					call rho_exp(rho(:,:,2),rho(:,:,1),tt)
				end if

				evops(1,1)%Uee(:,:,a,b,t_index/gt(2)+1) = rho(:,:,2)
			end if

		end do
		end if

		end do
		end do

		if(.not. evops_runge_kutta) then
			call clean_rho_exp()
		end if

		DEALLOCATE(rho)
		DEALLOCATE(rho0)
		DEALLOCATE(drho)



	    do a=1, size(evops(1,1)%UfeS,1)
	    do b=1, size(evops(1,1)%UfeS,2)
	    !do c=1, size(evops(1,1)%Ufe,3)
	    !do d=1, size(evops(1,1)%Ufe,4)
	    c = a
	    d = b

	    write(*,*) '-----------', a,b,en2(a), en(b), en2(a) - en(b) - rwa

	    do t=1, size(evops(1,1)%UfeS,3)

	    	evops(1,1)%UfeS(a,b,t) = exp(-(0,1)*(t-1)*dt*gt(1) &
	    		* (en2(a) - en(b) - rwa) )

	    	!evops(1,1)%UfeS(a,b,t) = evops(1,1)%UfeS(a,b,t) * exp(-gg_tr1((t-1)*gt(1)+1))
	    	!evops(1,1)%UfeS(a,b,t) = 1.0_dp
	    	!cycle

		!if(t > 1) then
	    	if(a <= Ncar) then
	    		if(b <= Ncar) then
	    			evops(1,1)%UfeS(a,b,t) = evops(1,1)%UfeS(a,b,t) * 0
	    		elseif(b <= 2*Ncar) then
	    			evops(1,1)%UfeS(a,b,t) = evops(1,1)%UfeS(a,b,t) * exp(-gg_tr3((t-1)*gt(1)+1))
	    		elseif(b <= 3*Ncar) then
	    			evops(1,1)%UfeS(a,b,t) = evops(1,1)%UfeS(a,b,t) * 0
	    		elseif(b <= 4*Ncar) then
	    			evops(1,1)%UfeS(a,b,t) = evops(1,1)%UfeS(a,b,t) * 0
	    		end if

	    	elseif(a <= 2*Ncar) then
	    		if(b <= Ncar) then
	    			evops(1,1)%UfeS(a,b,t) = evops(1,1)%UfeS(a,b,t) * 0
	    		elseif(b <= 2*Ncar) then
	    			evops(1,1)%UfeS(a,b,t) = evops(1,1)%UfeS(a,b,t) * exp(-gg_tr4((t-1)*gt(1)+1))
	    		elseif(b <= 3*Ncar) then
	    			evops(1,1)%UfeS(a,b,t) = evops(1,1)%UfeS(a,b,t) * exp(-gg_tr2((t-1)*gt(1)+1))
	    		elseif(b <= 4*Ncar) then
	    			evops(1,1)%UfeS(a,b,t) = evops(1,1)%UfeS(a,b,t) * 0
	    		end if
	    	elseif(a <= 3*Ncar) then
	    		if(b <= Ncar) then
	    			evops(1,1)%UfeS(a,b,t) = evops(1,1)%UfeS(a,b,t) * exp(-gg_tr4((t-1)*gt(1)+1))
	    		elseif(b <= 2*Ncar) then
	    			evops(1,1)%UfeS(a,b,t) = evops(1,1)%UfeS(a,b,t) * 0
	    		elseif(b <= 3*Ncar) then
	    			evops(1,1)%UfeS(a,b,t) = evops(1,1)%UfeS(a,b,t) * exp(-gg_tr1((t-1)*gt(1)+1))
	    		elseif(b <= 4*Ncar) then
	    			evops(1,1)%UfeS(a,b,t) = evops(1,1)%UfeS(a,b,t) * 0
	    		end if
	    	elseif(a <= 4*Ncar) then
	    		if(b <= Ncar) then
	    			evops(1,1)%UfeS(a,b,t) = evops(1,1)%UfeS(a,b,t) * exp(-gg_tr6((t-1)*gt(1)+1))
	    		elseif(b <= 2*Ncar) then
	    			evops(1,1)%UfeS(a,b,t) = evops(1,1)%UfeS(a,b,t) * 0
	    		elseif(b <= 3*Ncar) then
	    			evops(1,1)%UfeS(a,b,t) = evops(1,1)%UfeS(a,b,t) * 0
	    		elseif(b <= 4*Ncar) then
	    			evops(1,1)%UfeS(a,b,t) = evops(1,1)%UfeS(a,b,t) * exp(-gg_tr1((t-1)*gt(1)+1))
	    		end if
	    	elseif(a <= 5*Ncar) then
	    		if(b <= Ncar) then
	    			evops(1,1)%UfeS(a,b,t) = evops(1,1)%UfeS(a,b,t) * exp(-gg_tr5((t-1)*gt(1)+1))
	    		elseif(b <= 2*Ncar) then
	    			evops(1,1)%UfeS(a,b,t) = evops(1,1)%UfeS(a,b,t) * 0
	    		elseif(b <= 3*Ncar) then
	    			evops(1,1)%UfeS(a,b,t) = evops(1,1)%UfeS(a,b,t) * 0
	    		elseif(b <= 4*Ncar) then
	    			evops(1,1)%UfeS(a,b,t) = evops(1,1)%UfeS(a,b,t) * 0
	    		end if
	    	end if
	    !end if

	    evops(1,1)%UfeS(a,b,t) = conjg(evops(1,1)%UfeS(a,b,t))

	    !end do
	    !end do
	    end do
	    end do
	    end do

	    DEALLOCATE(iblocks(1,1)%sblock%dd)
	    DEALLOCATE(iblocks(1,1)%sblock%dx)
	    DEALLOCATE(iblocks(1,1)%sblock%dy)
	    DEALLOCATE(iblocks(1,1)%sblock%dz)
	    DEALLOCATE(iblocks(1,1)%sblock%en)
	    DEALLOCATE(iblocks(1,1)%sblock%J)

	    DEALLOCATE(iblocks(1,1)%eblock%dd)
	    DEALLOCATE(iblocks(1,1)%eblock%dx)
	    DEALLOCATE(iblocks(1,1)%eblock%dy)
	    DEALLOCATE(iblocks(1,1)%eblock%dz)
	    DEALLOCATE(iblocks(1,1)%eblock%en)
	    DEALLOCATE(iblocks(1,1)%eblock%eng)
	    DEALLOCATE(iblocks(1,1)%eblock%en_2)
	    DEALLOCATE(iblocks(1,1)%eblock%SS)
	    DEALLOCATE(iblocks(1,1)%eblock%S1)

	    DEALLOCATE(iblocks(1,1)%eblock%dd_2)
	    DEALLOCATE(iblocks(1,1)%eblock%dx_2)
	    DEALLOCATE(iblocks(1,1)%eblock%dy_2)
	    DEALLOCATE(iblocks(1,1)%eblock%dz_2)
	    DEALLOCATE(iblocks(1,1)%eblock%SS_2)
	    DEALLOCATE(iblocks(1,1)%eblock%S1_2)


	    ALLOCATE(iblocks(1,1)%sblock%dd,(Nvve,Nvvg))
	    ALLOCATE(iblocks(1,1)%sblock%dx,(Nvve,Nvvg))
	    ALLOCATE(iblocks(1,1)%sblock%dy,(Nvve,Nvvg))
	    ALLOCATE(iblocks(1,1)%sblock%dz,(Nvve,Nvvg))
	    ALLOCATE(iblocks(1,1)%sblock%en,(Nvve))
	    ALLOCATE(iblocks(1,1)%sblock%J,(Nvve,Nvve))

	    ALLOCATE(iblocks(1,1)%eblock%dd,(Nvve,Nvvg))
	    ALLOCATE(iblocks(1,1)%eblock%dx,(Nvve,Nvvg))
	    ALLOCATE(iblocks(1,1)%eblock%dy,(Nvve,Nvvg))
	    ALLOCATE(iblocks(1,1)%eblock%dz,(Nvve,Nvvg))
	    ALLOCATE(iblocks(1,1)%eblock%eng,(Nvvg))
	    ALLOCATE(iblocks(1,1)%eblock%en,(Nvve))
	    ALLOCATE(iblocks(1,1)%eblock%en_2,(Nvvf))
	    ALLOCATE(iblocks(1,1)%eblock%SS,(Nvve,Nvve))
	    ALLOCATE(iblocks(1,1)%eblock%S1,(Nvve,Nvve))

	    ALLOCATE(iblocks(1,1)%eblock%dd_2,(Nvvf,Nvve))
	    ALLOCATE(iblocks(1,1)%eblock%dx_2,(Nvvf,Nvve))
	    ALLOCATE(iblocks(1,1)%eblock%dy_2,(Nvvf,Nvve))
	    ALLOCATE(iblocks(1,1)%eblock%dz_2,(Nvvf,Nvve))
	    ALLOCATE(iblocks(1,1)%eblock%SS_2,(Nvvf,Nvvf))
	    ALLOCATE(iblocks(1,1)%eblock%S1_2,(Nvvf,Nvvf))

	    iblocks(1,1)%sblock%dd = dd
	    iblocks(1,1)%sblock%dx = dx
	    iblocks(1,1)%sblock%dy = dy
	    iblocks(1,1)%sblock%dz = dz
	    iblocks(1,1)%sblock%en = en
	    iblocks(1,1)%sblock%J = J

		! transform them here
		!! transformation is identical, we do not care for excitons

	    iblocks(1,1)%eblock%dd = dd
	    iblocks(1,1)%eblock%dx = dx
	    iblocks(1,1)%eblock%dy = dy
	    iblocks(1,1)%eblock%dz = dz
	    iblocks(1,1)%eblock%eng = eng
	    iblocks(1,1)%eblock%en = en
	    iblocks(1,1)%eblock%en_2 = en2
	    iblocks(1,1)%eblock%SS = 0.0_dp
	    iblocks(1,1)%eblock%S1 = 0.0_dp

	    iblocks(1,1)%eblock%dd_2 = dd2
	    iblocks(1,1)%eblock%dx_2 = dx2
	    iblocks(1,1)%eblock%dy_2 = dy2
	    iblocks(1,1)%eblock%dz_2 = dz2
	    iblocks(1,1)%eblock%SS_2 = 0.0_dp
	    iblocks(1,1)%eblock%S1_2 = 0.0_dp

	    do i=1, size(iblocks(1,1)%eblock%SS,1)
	    	iblocks(1,1)%eblock%SS(i,i) = 1.0_dp
	    	iblocks(1,1)%eblock%S1(i,i) = 1.0_dp
	    end do

	    do i=1, size(iblocks(1,1)%eblock%SS_2,1)
	    	iblocks(1,1)%eblock%SS_2(i,i) = 1.0_dp
	    	iblocks(1,1)%eblock%S1_2(i,i) = 1.0_dp
	    end do


		DEALLOCATE(dd)
		DEALLOCATE(dd2)
		DEALLOCATE(dx)
		DEALLOCATE(dx2)
		DEALLOCATE(dy)
		DEALLOCATE(dy2)
		DEALLOCATE(dz)
		DEALLOCATE(dz2)
		DEALLOCATE(J)
		DEALLOCATE(eng)
		DEALLOCATE(en)
		DEALLOCATE(en2)


	end subroutine prepare_very_specific_site_ops

	subroutine prepare_rho_exp(rho0)
		complex(dpc), dimension(:,:), intent(in) 	:: rho0
		complex(dpc), dimension(size(rho0,1),size(rho0,2)) 	:: xxx, rho
		real(dp) :: t
		integer(i4b) :: i,j,k,l, first_N, second_N

		first_N = size(rho0,1)
		second_N = size(rho0,2)

		ALLOCATE(rhoexpXX,(first_N*second_N,first_N*second_N))

		do k=1, first_N
		do l=1, second_N

		xxx = 0.0_dp
		xxx(k,l) = 1.0_dp

		call derivative_rates(t, xxx, rho)

		do i=1, first_N
		do j=1, second_N

		rhoexpXX(SUPERINDEX_FROM_K_L(i,j,second_N),SUPERINDEX_FROM_K_L(k,l,second_N)) = &
			rho(i,j)

		end do
		end do

		end do
		end do

	end subroutine

	subroutine clean_rho_exp()
		DEALLOCATE(rhoexpXX)
	end subroutine

	subroutine rho_exp(rho,rho0,t)
		complex(dpc), dimension(:,:), intent(in) 	:: rho0
		complex(dpc), dimension(:,:), intent(out) 	:: rho
		complex(dpc), dimension(:,:), allocatable 	:: XX, xxx
		complex(dpc), dimension(:), allocatable		:: YY
		real(dp), intent(in) :: t
		integer(i4b) :: i,j,k,l, first_N, second_N

		first_N = size(rho0,1)
		second_N = size(rho0,2)

		ALLOCATE(XX,(first_N*second_N,first_N*second_N))
		ALLOCATE(xxx,(first_N,second_N))
		ALLOCATE(YY,(first_N*second_N))

		do k=1, first_N
		do l=1, second_N

		YY(SUPERINDEX_FROM_K_L(k,l,second_N)) = rho0(k,l)

		end do
		end do

		XX = rhoexpXX * t

		call matrix_exp(XX)

		YY = matmul(XX,YY)

		do k=1, first_N
		do l=1, second_N

		rho(k,l) = YY(SUPERINDEX_FROM_K_L(k,l,second_N))

		end do
		end do

		DEALLOCATE(XX)
		DEALLOCATE(xxx)
		DEALLOCATE(YY)
	end subroutine rho_exp

	subroutine derivative_rates(t, y, dydt)
 		real(dp), intent(in) :: t
 		complex(dpc), dimension(:,:), intent(in) 	:: y
 		complex(dpc), dimension(:,:), intent(out) 	:: dydt

 		call derivative_rates_inner(t, y, dydt)

 	end subroutine derivative_rates

 	subroutine read_derivative_rates_from_file()
 		character(len=512) :: buff, buff2
 		real(dp) :: num
 		integer :: n

	 	open(unit=11,file=trim(trim(out_dir)//trim('/rates.dat')) , status='old', action='read', err=42)

 		k_et2 = -1.0_dp
 		k_et1 = -1.0_dp
 		k_ic21 = -1.0_dp
 		k_rate = -1.0_dp
 		k_sink = -1.0_dp


		do n=1,5
	 		read(11,*) buff, num
	 		write(buff2,*) 'reading derivative rates: ', trim(buff), num
	 		call print_log_message(trim(buff2),5)

		 	if(trim(buff) == 'k_et1') then
		 		k_et1 = num
		 	elseif(trim(buff) == 'k_et2') then
		 		k_et2 = num
		 	elseif(trim(buff) == 'k_ic21') then
		 		k_ic21 = num
		 	elseif(trim(buff) == 'k_sink') then
		 		k_sink = num
		 	elseif(trim(buff) == 'k_rate') then
		 		k_rate = num
		 	else
		 		write(*,*) 'ATTEMPT READING DERIVATIVE RATES, SOMETHING STRANGE HERE!: ', trim(buff), num
		 		stop
	 		end if
	 	end do

	 	if(k_et2 < 0 .or. k_et1 < 0 .or. k_ic21 < 0 .or. k_rate < 0 .or. k_sink < 0) then
	 		write(*,*) 'ATTEMPT READING DERIVATIVE RATES, INCOMPLETE'
		 	stop
	 	end if

		return

42		call print_error_message(-1,"out/rates.dat not found or incorrect format")
 	end subroutine

	subroutine derivative_rates_inner(t, y, dydt)
 		real(dp), intent(in) :: t
 		complex(dpc), dimension(:,:), intent(in) 	:: y
 		complex(dpc), dimension(:,:), intent(out) 	:: dydt

 		integer(i4b) :: a,b, Ncar, aa,bb, c,d, cc,dd

 		dydt = 0.0_dp

! 		k_et1 = 1/200.0_dp
!  		k_rate = 0.0_dp

! 		k_et2 = 0.0_dp+ 1/50.0_dp
! 		k_ic21 = 1/300.0_dp
! 		k_sink = 1/40.0_dp
! 		stop

 		!k_et2 = 0.0_dp
 		!k_et1 = 0.0_dp
 		!k_ic21 = 0.0_dp
 		!k_rate = 0.0_dp
 		!k_sink = 0.0_dp

 		!Ncar = current_s_block%QHO_lvls(1)
		Ncar = 6 ! LHQO

		! all coherences die at constant rate + coherent evolution
 		do a=1,size(y,1)
  		do b=1,size(y,1)
			dydt(a,b) = -k_rate*y(a,b) - cmplx(0,1)*(en(a)-en(b))*y(a,b)
		end do

			dydt(a,a) = 0.0_dp ! populations do not evolve
		end do

		! Kic12
		!do a=Ncar+1,2*Ncar
 		do c=1,Ncar
 		do d=1,Ncar

		! in this model, coherences are not transfered and vibrational populations
		! remain unchanged, only slowly transfer to S1
		if(c == d) then 	! population transfer with drho_jj/dt = - K_ij rho_ii

		dydt(c+Ncar,c+Ncar) = dydt(c+Ncar,c+Ncar) + k_ic21*y(c,c) !/Ncar
		dydt(c,c) = dydt(c,c) - k_ic21*y(c,c) !/Ncar

		else				! coherences are lowered as drho_ik/dt = - (K_ij + K_kj) rho_ik / 2

		dydt(c,d) = dydt(c,d) - k_ic21*y(c,d) !/Ncar
		! coherences simply do not appear in S1

		end if

 		!end do
 		end do
 		end do

		! Ket2 -- S2 to Qx
		do a=2*Ncar+1,3*Ncar
		do b=2*Ncar+1,3*Ncar
 		do c=1,Ncar
 		do d=1,Ncar

 		aa = a-2*Ncar
 		bb = b-2*Ncar
 		cc = c
 		dd = d

		! coherences survive but should not have effect, because of trace
		dydt(a,b) = dydt(a,b) + K_et2*LQHO(aa,"DD",cc)*LQHO(bb,"DD",dd)*y(c,d) !/Ncar
		dydt(c,d) = dydt(c,d) - K_et2*LQHO(aa,"DD",cc)*LQHO(bb,"DD",dd)*y(c,d) !/Ncar

 		end do
 		end do
 		end do
 		end do

		! Ket1 -- S1 to Qy
		do a=3*Ncar+1,4*Ncar
		do b=3*Ncar+1,4*Ncar
 		do c=1+Ncar,2*Ncar
 		do d=1+Ncar,2*Ncar

 		aa = a-3*Ncar
 		bb = b-3*Ncar
 		cc = c-Ncar
 		dd = d-Ncar

		! coherences survive but should not have effect, because of trace
		dydt(a,b) = dydt(a,b) + K_et1*LQHO(aa,"DD",cc)*LQHO(bb,"DD",dd)*y(c,d) !/Ncar
		dydt(c,d) = dydt(c,d) - K_et1*LQHO(aa,"DD",cc)*LQHO(bb,"DD",dd)*y(c,d) !/Ncar

 		end do
 		end do
 		end do
 		end do

 		!Ksink -- Qx to Qy
 		do a=2*Ncar+1,3*Ncar
 		do b=2*Ncar+1,3*Ncar

 		aa = a - 2*Ncar
 		bb = b - 2*Ncar

		dydt(a,b) = dydt(a,b)                     - K_sink*y(a,b) !/Ncar
		dydt(a+Ncar,b+Ncar) = dydt(a+Ncar,b+Ncar) + K_sink*y(a,b) !/Ncar

 		end do
 		end do


 	end subroutine derivative_rates_inner


	real(dp) function rho0_gg_element(ig,jg) result(res)
		integer(i4b), intent(in) :: ig,jg
		real(dp) :: ZZ
		integer(i4b) :: i

		if(ig /= jg) then
			res = 0.0_dp
			return
		end if

		ZZ = 0.0_dp
		do i=1, size(iblocks(1,1)%eblock%eng)
			ZZ = ZZ + exp(-iblocks(1,1)%eblock%eng(i) / temp / kB_intK)
		end do

		res = exp(-iblocks(1,1)%eblock%eng(ig) / temp / kB_intK)/ZZ
	end function rho0_gg_element

	real(dp) function pulse_spectrum(energy) result(spect)
		real(dp), intent(in) :: energy
		real(dp), dimension(:), allocatable :: tmp, tmp2
		real(dp) :: E_strength, linear, energy_inv_cm
		integer(i4b) :: n, i

		energy_inv_cm = Energy_internal_to_cm * energy
		spect = 0.0_dp

		open(unit=32,file=trim(trim(out_dir)//trim('/LO.dat')) , status='old', action='read', err=42)
    	read(32,*) n, E_strength

    	write(*,*) 'NACTENO',n, E_strength

     	ALLOCATE(tmp,(n))
    	ALLOCATE(tmp2,(n))

    	read(32,*) tmp
    	read(32,*) tmp2

    	tmp2 = tmp2 / maxval(tmp2)

    	close(32)

    	write(*,*) 'eeenergy',energy_inv_cm

    	do i=1, n-1
    		if(tmp(i+1) > energy_inv_cm .and. tmp(i) <= energy_inv_cm) then
    			linear = (energy_inv_cm-tmp(i))/(tmp(i+1)-tmp(i))
    			spect = sqrt((linear*(tmp2(i+1)-tmp2(i))+tmp2(i)))
    			    	write(*,*) 'eeenergy',energy_inv_cm, tmp(i), spect
    			exit
    		end if
    	end do

     	DEALLOCATE(tmp)
    	DEALLOCATE(tmp2)

		return

42		write(*,*) 'LO.dat not found'
		spect = 1.0_dp
	end function pulse_spectrum


 end module nakajima_zwanzig_shared
