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

        ! non-rephasing signal
        call fft_row_2D(SNR0,1,1)
        call fft_shift_2D(SNR0)

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


    subroutine save_2D_tot(data,filename)
        real, dimension(:,:), intent(in) :: data
        character(len=*), intent(in) :: filename
        character(len=20)	:: buffer, buffer2, buffer3 ! Origin-friendly output
        integer				:: log_max_data_value, f1, f2

       ! real(kind=1.0e

        write(buffer,'(i4)') Nt1

        log_max_data_value = int(log(1e-10+maxval(abs(data)))/log(10.0)) + 1
        f2 = min(max(8 - log_max_data_value,0),9)
        f1 = min(max(log_max_data_value + 4 + f2,13),99)
        write(buffer2,'(i2)') f1
        write(buffer3,'(i1)') f2

        buffer = '('//trim(buffer)//'f'//trim(buffer2)//'.'//trim(buffer3)//')'

        open(unit=11,file=filename,form="unformatted")
!        do i = 1, Nt1
!            do j = 1, Nt1

!               write(11,buffer) data(i,:) !real(SR0(i,:)+SNR0(i,:))
               write(unit=11) data

!            end do
!        end do
        close(11)

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

    subroutine save_2D_limits(dom,NFFT,filename)
        real(dp), intent(in) :: dom
        integer, intent(in) :: NFFT
        character(len=*), intent(in) :: filename
        real(dp) :: llim, ulim
        real(dp) :: width

        width = (real(NFFT,dp)/2.0_dp)*dom

        llim = rwa - width
        ulim = rwa + width

        open(unit=11,file=filename)
        write(11,*) llim*Energy_internal_to_cm, ulim*Energy_internal_to_cm, rwa*Energy_internal_to_cm
        close(11)

    end subroutine


    subroutine clean_twod

        call clean_signal3()

    end subroutine clean_twod

end module twod
