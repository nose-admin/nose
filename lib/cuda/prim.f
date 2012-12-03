# 1 "prim.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "prim.F90"
program testprim







   use iso_c_binding
   implicit none


    ! iso_c_binding interface
    interface
  subroutine primitive_GPU(Nx,Nw,x, y, w, d, pr, pi)
   use iso_c_binding
   bind(c,name="primitive_GPU_")
   type(c_ptr), value :: Nx,Nw
   type(c_ptr), value :: x,y,w,d,pr,pi
  end subroutine
  subroutine test_int(n)
   use iso_c_binding
   bind(c,name="test_int_")
   type(c_ptr), value :: n
  end subroutine
  subroutine test_float(x)
   use iso_c_binding
   bind(c,name="test_float_")
   type(c_ptr), value :: x
  end subroutine
 end interface
# 43 "prim.F90"
! external primitive_GPU


   integer(C_INT),target :: Nx
   integer(C_INT),target :: Nw
   real(C_FLOAT), dimension(:), allocatable,target :: x, y, d, pr, pi
   real(C_FLOAT), dimension(:), allocatable,target :: w






   integer :: i

   Nx = 1000
   Nw = 42

   allocate(x(Nx),y(Nx),d(Nx),pr(Nx*Nw),pi(Nx*Nw))
   allocate(w(Nw))

   do i = 1, Nx
      x(i) = (i-1)*0.1
      y(i) = sin(x(i))
   end do

   do i = 1, Nw
      w(i) = 10.0d0
   end do

   d = 0.0d0
   pr = 0.0d0
   pi = pr

   d(1) = 2.0
   d(2) = 1.0

   print *, "Calling C ..."

   print *, Nx
   call test_int(c_loc(Nx))
   print *, Nx
   call test_float(c_loc(x(2)))
   call primitive_GPU(c_loc(Nx),c_loc(Nw),c_loc(x),c_loc(y),&
                       c_loc(w),c_loc(d),c_loc(pr),c_loc(pi))
# 96 "prim.F90"
   print *, pr(1), pr(2)
   print *, pi(1), pi(2)



   print *, "...finished"
# 110 "prim.F90"
end program testprim
