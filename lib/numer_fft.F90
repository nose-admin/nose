!
! Fourier transform module
!
!
!
module numer_fft

    use std_types

    implicit none


    interface fft_row
    	! sum r=1,N exp(2 pi i isign (r-1)(s-1) / n) * u_r ;
    	! for inversion divide by N and use opposite isign
    	! w_max = 2 pi N / Tmax for it
        subroutine fft_row_sp(data,isign)
            use std_types
            complex(spc), dimension(:,:), intent(inout) :: data
            integer(i4b), intent(in)                    :: isign
        end subroutine fft_row_sp
        subroutine fft_row_dp(data,isign)
            use std_types
            complex(dpc), dimension(:,:), intent(inout) :: data
            integer(i4b), intent(in)                    :: isign
        end subroutine fft_row_dp
        module procedure fft_row_matrix_dp
    end interface

    interface fft_shift_row
        module procedure fft_shift_row_dpc
    end interface

    interface fft_regularize_row
        module procedure fft_regularize_row_dpc
    end interface

    interface fft_row_regularized
    	module procedure fft_row_regularized_dpc, fft_row_regularized_matrix_dpc
    end interface

    interface fft_row_convolution
    	module procedure fft_row_convolution_dpc
    end interface

	private::fft_regularize_row_dpc
	private::fft_regularize_row
	private::fft_row_regularized_dpc
	private::fft_row_regularized_matrix_dpc
	private::fft_shift_row_dpc
	private::fft_row_sp
	private::fft_row_dp
	private::fft_row_matrix_dp
	private::fft_row_2D_dpc
	private::fft_row_2D2_dpc
	private::fft_shift_2D_dpc
	private::fft_row_convolution_dpc



    !
    !
    ! 2D Fourier transform
    !
    ! Literature:
    !
    interface fft_row_2D
        module procedure fft_row_2D_dpc, fft_row_2D2_dpc
    end interface

    interface fft_shift_2D
        module procedure fft_shift_2D_dpc
    end interface

contains

    !
    ! Two-dimensional Fourier transform
    !
    subroutine fft_row_2D_dpc(data,isign)
        complex(DPC), dimension(:,:), intent(INOUT) :: data
        integer(I4B), intent(IN) :: isign
        complex(DPC), dimension(size(data,2),size(data,1)) :: temp
        call fft_row(data,isign)
        temp=transpose(data)
        call fft_row(temp,isign)
        data=transpose(temp)
    end subroutine fft_row_2D_dpc

    !
    ! Two-dimensional Fourier transform
    !
    subroutine fft_row_2D2_dpc(data,isign1,isign2)
        complex(DPC), dimension(:,:), intent(INOUT) :: data
        integer(I4B), intent(IN) :: isign1, isign2
        complex(DPC), dimension(size(data,2),size(data,1)) :: temp
        ! first the FFT is done in the second index
        call fft_row(data,isign2)
        temp=transpose(data)
        ! then in the first
        call fft_row(temp,isign1)
        data=transpose(temp)
    end subroutine fft_row_2D2_dpc

	!
	! FFT-shift in both indexes
	!
    subroutine fft_shift_2D_dpc(data)
        complex(DPC), dimension(:,:), intent(INOUT) :: data
        ! local
        complex(dpc), dimension(1:size(data,1),1:size(data,2)) :: dt
        integer :: i,j, N1, N2

        N1 = size(data,1)
        N2 = size(data,2)

        if(mod(N1,2) > 0 .or. mod(N2,2) > 0) then
        	write(*,*) 'error in fft_shift_2D_dpc - argument not 2-divisible'
        	stop
        endif

        dt(N1/2+1:N1,N2/2+1:N2) = data(1:N1/2,1:N2/2)
        dt(1:N1/2,1:N2/2)       = data(N1/2+1:N1,N2/2+1:N2)
        dt(1:N1/2,N2/2+1:N2)    = data(N1/2+1:N1,1:N2/2)
        dt(N1/2+1:N1,1:N2/2)    = data(1:N1/2,N2/2+1:N2)

        data = dt

    end subroutine fft_shift_2D_dpc

	!
	! FFT-shift in row index
	!
    subroutine fft_shift_row_dpc(data)
        complex(DPC), dimension(:,:), intent(INOUT) :: data
        ! local
        complex(dpc), dimension(1:size(data,1),1:size(data,2)) :: dt
        integer :: i, j, N2

        N2 = size(data,2)

        if(mod(N2,2) > 0) then
        	write(*,*) 'error in fft_shift_row_dpc - argument not 2-divisible'
        	stop
        endif

        dt(:,N2/2+1:N2) = data(:,1:N2/2)
        dt(:,1:N2/2)    = data(:,N2/2+1:N2)

        data = dt

    end subroutine fft_shift_row_dpc

    ! regularized FFTs

	subroutine fft_regularize_row_dpc(data, half_time_in_Tmax, isign)
		! Tmax = size(data,2)/2*dt

		double complex, dimension(:,:), intent(inout) 	:: data
		double precision, intent(in) 						:: half_time_in_Tmax
		double precision 	:: Tmax_ratio
		integer, intent(in)								:: isign
		integer 			:: i, j
		integer				:: N1, N2

		N1 = size(data,1)
		N2 = size(data,2)

        if(mod(N2,2) > 0) then
        	write(*,*) 'error in fft_regularize_row_dpc - argument not 2-divisible'
        	stop
        endif

		do i=1, N1
		do j=1, N2
			if(j <= N2/2) then
				Tmax_ratio = ((REAL(j)-1)/(N2/2))
			else
				Tmax_ratio = ((REAL(j)-1-N2)/(N2/2))
			endif

			data(i,j) = data(i,j)*exp(isign/half_time_in_Tmax*Tmax_ratio)
		end do
		end do

	end subroutine fft_regularize_row_dpc

	function half_time_in_Tmax_to_regularization_alpha(half_time_in_Tmax,Tmax) result(regularization_alpha)
		double precision, intent(in)	:: half_time_in_Tmax
		double precision, intent(in)	:: Tmax
		double precision				:: regularization_alpha

		regularization_alpha = 1/(half_time_in_Tmax*Tmax)
	end function half_time_in_Tmax_to_regularization_alpha

	subroutine fft_row_regularized_dpc(data,isign,half_time_in_Tmax)
		complex(dpc), dimension(:,:), intent(inout) 	:: data
		integer(i4b), intent(in)                    	:: isign
		double precision, intent(in)					:: half_time_in_Tmax

		if(isign == 1) then
			call fft_regularize_row(data, half_time_in_Tmax, -isign)
			call fft_row(data, isign)
		elseif(isign == -1) then
			call fft_row(data, isign)
			call fft_regularize_row(data, half_time_in_Tmax, -isign)
			data = data/size(data,2)
		else
			write(*,*) 'error in fft_row_regularized_dpc'
			stop
		end if

	end subroutine fft_row_regularized_dpc

	subroutine fft_row_regularized_matrix_dpc(data,isign,half_time_in_Tmax)
		complex(dpc), dimension(:,:,:), intent(inout) 	:: data
		integer(i4b), intent(in)  		                  	:: isign
		double precision, intent(in)						:: half_time_in_Tmax
		integer(i4b) :: i

		do i=1,size(data,1)
			call fft_row_regularized_dpc(data(i,:,:),isign,half_time_in_Tmax)
		end do

	end subroutine fft_row_regularized_matrix_dpc

	! convolution

	subroutine fft_row_convolution_dpc(to, f1, f2, half_time_in_Tmax)
		complex(dpc), dimension(:,:), intent(in) 		:: f1, f2
		complex(dpc), dimension(:,:), intent(out)		:: to
		complex(dpc), dimension(size(to,1),size(to,2)):: tmp
		double precision, intent(in)					:: half_time_in_Tmax

		if(.not.(size(f1,1)==size(f2,1) .and. size(f1,2)==size(f2,2) &
			.and. size(f1,2)==size(to,2) .and. size(f1,2)==size(to,2))) then
			write(*,*) 'dimension error in fft_row_convolution_dpc()'
			stop
		end if

		to  = f1
		tmp = f2

		call fft_row_regularized(to,   1, half_time_in_Tmax)
		call fft_row_regularized(tmp,  1, half_time_in_Tmax)
		to = to*tmp
		call fft_row_regularized(to, -1, half_time_in_Tmax)

	end subroutine fft_row_convolution_dpc

	subroutine fft_row_matrix_dp(data,isign)
		complex(dpc), dimension(:,:,:), intent(inout) 	:: data
		integer(i4b), intent(in)                    		:: isign
		integer(i4b)										:: i

		do i=1,size(data,1)
			call fft_row(data(i,:,:),isign)
		end do

	end subroutine fft_row_matrix_dp

end module numer_fft
