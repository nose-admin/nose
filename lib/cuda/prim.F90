program testprim

#define NOSE_ISO_C_BINDING
#undef NOSE_ISO_C_BINDING

#ifdef USE_CUDA


#ifdef NOSE_ISO_C_BINDING
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
		subroutine test_ptr(x,n)
			use iso_c_binding
			bind(c,name="test_ptr_")
			type(c_ptr), value :: x
			type(c_ptr), value :: n
		end subroutine
	end interface
#else
   implicit none
   interface
      subroutine primitive_GPU(Nx,Nw,x, y, w, d, pr, pi)
      real*4, dimension(*) :: x, y, w, d, pi, pr
      integer              :: Nx, Nw
      end subroutine
   end interface
#endif


!   external primitive_GPU

#ifdef NOSE_ISO_C_BINDING
   integer(C_INT),target :: Nx
   integer(C_INT),target :: Nw
   real(C_FLOAT), dimension(:), allocatable,target  :: x, y, d, pr, pi
   real(C_FLOAT), dimension(:), allocatable,target  :: w
#else
   integer :: Nx
   integer :: Nw
   real*4, dimension(:), allocatable  :: x, y, d, pr, pi
   real*4, dimension(:), allocatable  :: w
#endif
   integer               :: i

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
#ifdef NOSE_ISO_C_BINDING
   print *, Nx
   call test_int(c_loc(Nx))
   print *, Nx
   call test_float(c_loc(x(2)))
   call primitive_GPU(c_loc(Nx),c_loc(Nw),c_loc(x),c_loc(y),&
                       c_loc(w),c_loc(d),c_loc(pr),c_loc(pi))

   Nx = 3
   call test_ptr(c_loc(x),c_loc(Nx))
   print *, "test ", x(3)

#else
   print *, Nx
   call test_int(Nx)
   print *, Nx
   call test_float(x(2))
   call primitive_GPU(Nx,Nw,x,y,w,d,pr,pi)

#endif
   print *, pr(1), pr(2)
   print *, pi(1), pi(2)



   print *, "...finished"
#endif

#ifndef USE_CUDA

    print *, "CUDA Parallelization is not used"

#endif

end program testprim
