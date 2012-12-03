!!
!! Ordinary Differential Equations Solvers
!!
module numer_ode

  implicit none

  !!
  !! 4-th order Runge-Kutta method
  !!
  interface ode_rk4

     subroutine rk4_sp(y,dydx,x,h,yout,derivs)
       use std_types
       real(SP), dimension(:), intent(IN) :: y,dydx
       real(SP), intent(IN) :: x,h
       real(SP), dimension(:), intent(OUT) :: yout
       interface
          subroutine derivs(x,y,dydx)
            use std_types
            real(SP), intent(IN) :: x
            real(SP), dimension(:), intent(IN) :: y
            real(SP), dimension(:), intent(OUT) :: dydx
          end subroutine derivs
       end interface
     end subroutine rk4_sp

     subroutine rk4_dp(y,dydx,x,h,yout,derivs)
       use std_types
       real(DP), dimension(:), intent(IN) :: y,dydx
       real(DP), intent(IN) :: x,h
       real(DP), dimension(:), intent(OUT) :: yout
       interface
          subroutine derivs(x,y,dydx)
            use std_types
            real(DP), intent(IN) :: x
            real(DP), dimension(:), intent(IN) :: y
            real(DP), dimension(:), intent(OUT) :: dydx
          end subroutine derivs
       end interface
     end subroutine rk4_dp
     
     subroutine rk4_dpc(y,dydx,x,h,yout,derivs)
       use std_types
       complex(DPC), dimension(:), intent(IN) :: y,dydx
       real(DP), intent(IN) :: x,h
       complex(DPC), dimension(:), intent(OUT) :: yout
       interface
          subroutine derivs(x,y,dydx)
            use std_types
            real(DP), intent(IN) :: x
            complex(DPC), dimension(:), intent(IN) :: y
            complex(DPC), dimension(:), intent(OUT) :: dydx
          end subroutine derivs
       end interface
     end subroutine rk4_dpc

     subroutine rk4_dpc_mat(y,dydx,x,h,yout,derivs)
       use std_types
       complex(DPC), dimension(:,:), intent(IN) :: y,dydx
       real(DP), intent(IN) :: x,h
       complex(DPC), dimension(:,:), intent(OUT) :: yout
       interface
          subroutine derivs(x,y,dydx) 
            use std_types
            real(DP), intent(IN) :: x
            complex(DPC), dimension(:,:), intent(IN) :: y
            complex(DPC), dimension(:,:), intent(OUT) :: dydx
          end subroutine derivs
       end interface
     end subroutine rk4_dpc_mat

     subroutine rk4_dpc_mat3(y,dydx,x,h,yout,derivs)
       use std_types
       complex(DPC), dimension(:,:,:), intent(IN)  :: y,dydx
       real(DP), intent(IN) :: x,h
       complex(DPC), dimension(:,:,:), intent(OUT) :: yout
       interface
          subroutine derivs(x,y,dydx) 
            use std_types
            real(DP), intent(IN) :: x
            complex(DPC), dimension(:,:,:), intent(IN)  :: y
            complex(DPC), dimension(:,:,:), intent(OUT) :: dydx
          end subroutine derivs
       end interface
     end subroutine rk4_dpc_mat3

  end interface

end module numer_ode

