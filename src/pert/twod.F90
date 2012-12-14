#include "util_allocation.h"

!
! Module containing subroutines for 2D spectroscopy
!
module twod

    use std_types
    use signal3
    use numer_fft

contains


    subroutine init_twod()

        integer :: i,j

        Nt1 = Nt(1)

        call init_signal3()

        ! rephasing signal
        call fft_row_2D(SR0,1,-1)
        call fft_shift_2D(SR0)

   		! rephasing signal
        call fft_row_2D(SSR2,1,-1)
        call fft_shift_2D(SSR2)

        ! non-rephasing signal
        call fft_row_2D(SNR0,1,1)
        call fft_shift_2D(SNR0)

        ! non-rephasing signal
        call fft_row_2D(SSR1,1,1)
        call fft_shift_2D(SSR1)

		! ESA
		call fft_row_2D(ESAR,1,-1)
		call fft_shift_2D(ESAR)
		call fft_row_2D(ESAN,1,1)
		call fft_shift_2D(ESAN)

        !stop

    end subroutine init_twod


!                open(unit=11,file=trim(file_join(out_dir,"pop_n_o_coh.dat")))
!
!                write(buffer,'(i3)') 2*((rcohs%N)**2) + 1
!                buffer = '('//trim(buffer)//'f12.8)'
!
!                do i = 1, Nt(1)
!
!                  write(11,buffer) (i-1)*gt(1)*dt, rcohs%RC(i,1:rcohs%N,1:rcohs%N)
!
!                end do
!
!               close(unit=11)


    subroutine save_2D_tot(data,filename,dom,NFFT,minOm,maxOm)
        real, dimension(:,:), intent(in) :: data
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: dom
        real(dp), optional, intent(in) :: minOm, maxOm

        character(len=20)	:: buffer, buffer2, buffer3 ! Origin-friendly output
        integer(i4b), intent(in) 	:: NFFT
        integer(i4b)				:: log_max_data_value, f1, f2, data_size, ii, jj
        real(dp)					:: width, llim, ulim, maxO, minO

        real(dp), dimension(size(data,1),size(data,2)) :: weigth,data2
        real(dp), allocatable, dimension(:,:) :: data3

        call load_LO_spectrum_weights(weigth,dom,NFFT)
		data2 = data*weigth

		! trim what we output according to frequencies
		if(present(minOm) .and. present(maxOm)) then
	        maxO = maxOm
    	    minO = minOm

			width = (real(NFFT,dp)/2.0_dp)*dom
			llim = (rwa - width)*Energy_internal_to_cm
			ulim = (rwa + width)*Energy_internal_to_cm

			if(maxO > ulim) then
				maxO = ulim
			end if
			if(minO < llim) then
				minO = llim
			end if

			data_size = int(size(data,1)*(maxO - minO)/(ulim-llim))

			if(data_size <= 0) then
				call print_error_message(-1, "chosen output frequencies create field of zero length")
				stop
			end if

			write(*,*) 'data size', data_size, maxO, minO, maxO - minO, (ulim-llim), size(data,1)

			ALLOCATE(data3,(data_size,data_size))

			data3 = 0.0_dp
			do i=1, size(data2,1)
			do j=1, size(data2,2)
				ii = int(i + (llim-minOm)*size(data2,1)/(ulim-llim) + 0.5)
				jj = int(j + (llim-minOm)*size(data2,1)/(ulim-llim) + 0.5)

				if(ii >= 1 .and. ii <= size(data3,1) .and. &
					jj >= 1 .and. jj <= size(data3,2)) then

					data3(ii,jj) = data2(i,j)
				end if
			end do
			end do
		else
			ALLOCATE(data3, (size(data,1),size(data,2)))
			data3 = data2
		end if

        write(buffer,'(i4)') Nt1

        log_max_data_value = int(log(1e-10+maxval(abs(data3)))/log(10.0)) + 1
        f2 = min(max(8 - log_max_data_value,0),9)
        f1 = min(max(log_max_data_value + 4 + f2,13),99)
        write(buffer2,'(i2)') f1
        write(buffer3,'(i1)') f2

        buffer = '('//trim(buffer)//'f'//trim(buffer2)//'.'//trim(buffer3)//')'

        open(unit=11,file=filename)
        !open(unit=11,file=filename,form="unformatted")
        do i = 1, size(data3,1)
!            do j = 1, Nt1

               write(11,buffer) data3(i,:)
!               write(unit=11) data

!            end do
        end do
        close(11)

        open(unit=12,file="twod_LO.dat")
        !open(unit=11,file=filename,form="unformatted")
        do i = 1, Nt1
!            do j = 1, Nt1

               write(12,buffer) weigth(i,:)
!               write(unit=11) data

!            end do
        end do
        close(12)

        DEALLOCATE(data3)

    end subroutine save_2D_tot

    subroutine save_2DRe_R(filename)
        character(len=*), intent(in) :: filename

        open(unit=11,file=filename)
        do i = 1, Nt1
        !    do j = 1, Nt1
               write(11,*) real(SR0(i,:))
        !    end do
        end do
        close(11)

    end subroutine save_2DRe_R

    subroutine save_2DRe_NR(filename)
        character(len=*), intent(in) :: filename

        open(unit=11,file=filename)
        do i = 1, Nt1
        !    do j = 1, Nt1
               write(11,*) real(SNR0(i,:))
        !    end do
        end do
        close(11)

    end subroutine save_2DRe_NR


    subroutine save_2DIm_R(filename)
        character(len=*), intent(in) :: filename

        open(unit=11,file=filename)
        do i = 1, Nt1
        !    do j = 1, Nt1
               write(11,*) aimag(SR0(i,:))
        !    end do
        end do
        close(11)

    end subroutine save_2DIm_R

    subroutine save_2DIm_NR(filename)
        character(len=*), intent(in) :: filename

        open(unit=11,file=filename)
        do i = 1, Nt1
        !    do j = 1, Nt1
               write(11,*) aimag(SNR0(i,:))
        !    end do
        end do
        close(11)

    end subroutine save_2DIm_NR

    subroutine save_2D_limits(dom,NFFT,filename,minOm,maxOm)
        real(dp), intent(in) :: dom
        integer, intent(in) :: NFFT
        character(len=*), intent(in) :: filename
        real(dp), intent(in), optional :: minOm, maxOm
        real(dp) :: llim, ulim
        real(dp) :: width

        width = (real(NFFT,dp)/2.0_dp)*dom

        llim = rwa - width
        ulim = rwa + width

        if(present(minOm) .and. present(maxOm)) then
        	llim = max(llim,minOm/Energy_internal_to_cm)
        	ulim = min(ulim,maxOm/Energy_internal_to_cm)
        end if

        write(*,*) '   ** -- ** limits',NFFT,width, llim, ulim, dom
        write(*,*) '   ** -- ** limits',NFFT,width*Energy_internal_to_cm, llim*Energy_internal_to_cm, ulim*Energy_internal_to_cm, dom


        open(unit=11,file=filename)
        write(11,*) llim*Energy_internal_to_cm, ulim*Energy_internal_to_cm, rwa*Energy_internal_to_cm
        close(11)

    end subroutine

    subroutine load_LO_spectrum_weights(ret,dom,NFFT)
    	real(dp), intent(out), dimension(:,:) :: ret
    	real(dp), intent(in) :: dom
        integer(i4b), intent(in) :: NFFT
!        character(len=*), intent(in) :: filename
        real(dp) :: llim, ulim, now_I, now_J
        real(dp) :: width
    	real(dp), allocatable, dimension(:) :: tmp, tmp2
    	integer(i4b) :: n, i, j, k
    	real(dp) :: linear, E_strength
    	character(len=25600) :: cbuff


        if(size(ret,1) /= size(ret,2)) then
        	call print_error_message(-1,"load_LO_spectrum_weights - non-square 2D spectrum given")
        end if

		!NFFT = size(ret,1)
        width = (real(NFFT,dp)/2.0_dp)*dom

        llim = rwa - width
        ulim = rwa + width

        write(*,*) '   ** -- ** ',NFFT,width, llim, ulim, dom
        write(*,*) '   ** -- ** ',NFFT,width*Energy_internal_to_cm, llim*Energy_internal_to_cm, ulim*Energy_internal_to_cm, dom


		write(cbuff,*) 'Trying LO data in file' ,trim(trim(out_dir)//trim('/LO.dat'))
		call print_log_message(trim(cbuff),5)
		call flush()

    	ret = 1.0_dp

    	open(unit=11,file=trim(trim(out_dir)//trim('/LO.dat')) , status='old', action='read', err=42)
    		read(11,*) n, E_strength

    		write(*,*) 'NACTENO',n, E_strength

 	    	ALLOCATE(tmp,(n))
    		ALLOCATE(tmp2,(n))

    		read(11,*) tmp
    		read(11,*) tmp2

    		write(cbuff,*) 'Laser strength read: ', tmp
    		call print_log_message(cbuff,5)
    		write(cbuff,*) 'Laser strength read: ', tmp2
    		call print_log_message(cbuff,5)

    		tmp = tmp/Energy_internal_to_cm
    		tmp2 = tmp2/maxval(tmp2)*E_strength

    		ret = 1.0_dp

    		do i=1,size(ret,1)
    		do j=1,size(ret,2)

    			now_I = REAL((i-1))/size(ret,1)*(ulim-llim)+llim
    			now_J = REAL((j-1))/size(ret,2)*(ulim-llim)+llim

    			!write(*,*) '    ***',i,j,now_I, now_J

    			do k=2,n
    				!write(*,*) '    *****',k
    				if(k == n) then
    				  ret(i,j) = 0.0_dp
    				  exit
    				end if
    				if(k == 2 .and. tmp(k) > now_I) then
    				  ret(i,j) = 0.0_dp
    				  exit
    				end if


    				if(tmp(k) > now_I) then
    					linear = (now_I-tmp(k-1))/(tmp(k)-tmp(k-1))
    					ret(i,j) = ret(i,j) * ( (linear*(tmp2(k)-tmp2(k-1))+tmp2(k-1)) )
    					exit
    				end if
    			end do

    			do k=2,n
    				if(k == n) then
    				  ret(i,j) = 0.0_dp
    				  exit
    				end if
    				if(k == 2 .and. tmp(k) > now_J) then
    				  ret(i,j) = 0.0_dp
    				  exit
    				end if

    				!write(*,*) '    *****',K
    				if(tmp(k) > now_J) then
    					linear = (now_J-tmp(k-1))/(tmp(k)-tmp(k-1))
    					ret(i,j) = ret(i,j) * sqrt(max(0.0_dp,linear*(tmp2(k)-tmp2(k-1))+tmp2(k-1)))
    					exit
    				end if
    			end do

    		end do
    		end do

	    	DEALLOCATE(tmp)
    		DEALLOCATE(tmp2)

    	close(11)

    	return


42		call print_log_message("out/LO.dat not found or incorrect format",5)

    end subroutine


    subroutine clean_twod

        call clean_signal3()

    end subroutine clean_twod

end module twod
