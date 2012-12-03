!
! Basic module for handling points in space
!
!
!
module util_space

  use std_types

  implicit none

  private

  !************************************************************************
  ! DECLARATIONS
  !************************************************************************
  real(dp), dimension(3,3)  :: RX,RY,RZ     ! rotation matrices

  integer(i4b) :: rotation_set


  public :: set_rotation
  public :: rotate_point
  public :: rotate_group

contains

  !
  ! Initionalization of the rotation matrices
  !
  subroutine set_rotation(thX,thY,thZ)
    real(dp), intent(in)  :: thX,thY,thZ

    RX = 0.0_dp
    RY = 0.0_dp
    RZ = 0.0_dp

    RX(1,1) = 1.0_dp
    RX(2,2) = cos(thX)
    RX(2,3) = -sin(thX)
    RX(3,2) = sin(thX)
    RX(3,3) = cos(thX)

    RY(1,1) = cos(thY)
    RY(1,3) = sin(thY)
    RY(2,2) = 1.0_dp
    RY(3,1) = -sin(thY)
    RY(3,3) = cos(thY)

    RZ(1,1) = cos(thZ)
    RZ(1,2) = -sin(thZ)
    RZ(2,1) = sin(thZ)
    RZ(2,2) = cos(thZ)
    RZ(3,3) = 1.0_dp

    rotation_set = 1

  end subroutine set_rotation

  subroutine rotate_along_X(dd)
    real(dp), dimension(:), intent(inout) :: dd
    
    if (size(dd) == 3) then
       dd = matmul(RX,dd)
    else
       write(*,*) "Error in: rotate_along_X; expected rank 3"
       stop
    end if

  end subroutine rotate_along_X

  subroutine rotate_along_Y(dd)
    real(dp), dimension(:), intent(inout) :: dd
    
    if (size(dd) == 3) then
       dd = matmul(RY,dd)
    else
       write(*,*) "Error in: rotate_along_X; expected rank 3"
       stop
    end if

  end subroutine rotate_along_Y

  subroutine rotate_along_Z(dd)
    real(dp), dimension(:), intent(inout) :: dd
    
    if (size(dd) == 3) then
       dd = matmul(RZ,dd)
    else
       write(*,*) "Error in: rotate_along_X; expected rank 3"
       stop
    end if

  end subroutine rotate_along_Z

  !
  ! Rotates along all three axes in specified order
  !
  subroutine rotate_point(ST,ND,RD,dd)
    character, intent(in) :: ST,ND,RD
    real(dp), dimension(:), intent(inout) :: dd

    character(len=1),dimension(3) :: ord
    integer(i4b) :: i

    if (rotation_set /= 1) then
       write(*,*) "Error: Rotation not initialized"
       stop 
    end if

    ord(1) = ST
    ord(2) = ND
    ord(3) = RD

    do i = 1, 3
       if (ord(i) == 'X') then
          call rotate_along_X(dd)
       else if (ord(i) == 'Y') then
          call rotate_along_Y(dd)
       else if (ord(i) == 'Z') then
          call rotate_along_Z(dd)
       else
          write(*,*) "Error: "
          stop
       end if
    end do

  end subroutine rotate_point
    
  !
  ! Rotates a group of points
  !
  subroutine rotate_group(ST,ND,RD,dd)
    character, intent(in) :: ST,ND,RD
    real(dp), dimension(:,:), intent(inout) :: dd

    character(len=1),dimension(3) :: ord
    integer(i4b) :: i, N, k

    if (rotation_set /= 1) then
       write(*,*) "Error: Rotation not initialized"
       stop 
    end if

    ord(1) = ST
    ord(2) = ND
    ord(3) = RD

    N = size(dd,2)
    
    do k = 1, N

       do i = 1, 3
          if (ord(i) == 'X') then
             call rotate_along_X(dd(:,k))
          else if (ord(i) == 'Y') then
             call rotate_along_Y(dd(:,k))
          else if (ord(i) == 'Z') then
             call rotate_along_Z(dd(:,k))
          else
             write(*,*) "Error: "
             stop
          end if
       end do
       
    end do

  end subroutine rotate_group


end module util_space
