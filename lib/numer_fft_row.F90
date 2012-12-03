subroutine fft_row_sp(data,isign)
  use std_types
!  use nrutil, only : assert,swap
  implicit none
  complex(SPC), dimension(:,:), intent(INOUT) :: data
  integer(I4B), intent(IN) :: isign
  integer(I4B) :: i,istep,j,m,mmax,n2
  integer :: n
  real(DP) :: theta
  complex(SPC), dimension(size(data,1)) :: temp
  complex(DPC) :: w,wp
  complex(SPC) :: ws
  n=size(data,2)
!  call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourrow_sp')
  if (.not.(iand(n,n-1)==0)) then
     write(*,*) "Error: n must be a power of 2 in fourrow_sp" 
     stop 
  end if
  n2=n/2
  j=n2
  do i=1,n-2
     if (j > i) then !call swap(data(:,j+1),data(:,i+1))
        temp = data(:,j+1)
        data(:,j+1) = data(:,i+1)
        data(:,i+1) = temp
     end if
     m=n2
     do
        if (m < 2 .or. j < m) exit
        j=j-m
        m=m/2
     end do
     j=j+m
  end do
  mmax=1
  do
     if (n <= mmax) exit
     istep=2*mmax
     theta=PI_D/(isign*mmax)
     wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
     w=cmplx(1.0_dp,0.0_dp,kind=dpc)
     do m=1,mmax
        ws=w
        do i=m,n,istep
           j=i+mmax
           temp=ws*data(:,j)
           data(:,j)=data(:,i)-temp
           data(:,i)=data(:,i)+temp
        end do
        w=w*wp+w
     end do
     mmax=istep
  end do
end subroutine fft_row_sp

subroutine fft_row_dp(data,isign)
  use std_types
!  use nrutil, only : assert,swap
  implicit none
  complex(DPC), dimension(:,:), intent(INOUT) :: data
  integer(I4B), intent(IN) :: isign
  integer(I4B) :: i,istep,j,m,mmax,n2
  integer :: n
  real(DP) :: theta
  complex(DPC), dimension(size(data,1)) :: temp
  complex(DPC) :: w,wp
  complex(DPC) :: ws
  n=size(data,2)
!  call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourrow_dp')
  if (.not.(iand(n,n-1)==0)) then
     write(*,*) "Error: n must be a power of 2 in fourrow_dp" 
     stop 
  end if
  n2=n/2
  j=n2
  do i=1,n-2
     if (j > i) then !call swap(data(:,j+1),data(:,i+1))
        temp = data(:,j+1)
        data(:,j+1) = data(:,i+1)
        data(:,i+1) = temp
     end if
     m=n2
     do
        if (m < 2 .or. j < m) exit
        j=j-m
        m=m/2
     end do
     j=j+m
  end do
  mmax=1
  do
     if (n <= mmax) exit
     istep=2*mmax
     theta=PI_D/(isign*mmax)
     wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
     w=cmplx(1.0_dp,0.0_dp,kind=dpc)
     do m=1,mmax
        ws=w
        do i=m,n,istep
           j=i+mmax
           temp=ws*data(:,j)
           data(:,j)=data(:,i)-temp
           data(:,i)=data(:,i)+temp
        end do
        w=w*wp+w
     end do
     mmax=istep
  end do
end subroutine fft_row_dp
