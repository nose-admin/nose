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

	implicit none

	public::N1_from_type
	public::N2_from_type
	public::operator_to_exc
	public::operator_from_exc
	public::superops_to_exc ! transforms the superoperator into exc picture directly

	interface superops_to_exc
		module procedure superops_to_exc_2indexed
		module procedure superops_to_exc_4indexed
	end interface

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
 			fromto = matmul(iblocks(1,1)%eblock%S1,fromto)
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
 			fromto = matmul(iblocks(1,1)%eblock%SS,fromto)
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


 end module nakajima_zwanzig_shared
