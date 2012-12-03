module numer_interp

	use numer_matrix ! propagates: std_types
	
	implicit none
	
	private :: get_abcd, tridag, locate
	
!	public :: spline, splint, primitive, intsplin

	interface spline
		module procedure spline_real
		module procedure spline_cmplx
		module procedure spline_matrix
	end interface
	interface splint
		module procedure splint_real
		module procedure splint_cmplx
	end interface

!	interface primitive
!		module procedure primitive_real
!		module procedure primitive_cmplx
!	end interface primitive


contains


    !
    ! Interpolation by spline between two points
    !
	function splint_real(xa,ya,y2a,x) result(ret)
		real(dp), dimension(:), intent(in) :: xa,ya,y2a
		real(dp), intent(in) :: x
		real(dp) :: ret
		
		integer(i4b) :: khi,klo,n
		real(dp) :: a,b,h
		
		n = size(xa)
		
		klo=max(min(locate(xa,x),n-1),1)
		khi=klo+1
		h=xa(khi)-xa(klo)
		!print *, klo,khi,h
		if (h == 0.0_dp) stop "Error in splint (1)"
		a = (xa(khi)-x)/h
		b = (x-xa(klo))/h
		ret = a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))* &
		 (h**2)/6.0_dp
	
	end function splint_real

	!
	! Cubic spline interpolation
	!
	subroutine spline_real(x,y,yp1,ypn,y2)
		real(dp), dimension(:), intent(in)  :: x,y
		complex(dpc), intent(in)                :: yp1, ypn
		real(dp), dimension(:), intent(out) :: y2 
	
		integer(i4b) :: n
		real(dp), dimension(size(x)) :: a,b,c,r
		n = size(x)
		c(1:n-1)=x(2:n)-x(1:n-1)
		r(1:n-1)=6.0_dp*((y(2:n)-y(1:n-1))/c(1:n-1))
		r(2:n-1)=r(2:n-1)-r(1:n-2)
		a(2:n-1)=c(1:n-2)
		b(2:n-1)=2.0_dp*(c(2:n-1)+a(2:n-1))
		b(1)=1.0_dp
		b(n)=1.0_dp
		
		if(abs(yp1) > 0.99e30_dp) then
			r(1) = 0.0_dp
			c(1) = 0.0_dp
		else
			r(1) = (3.0_dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
			c(1) = 0.5_dp
		end if
		
		if(abs(ypn) > 0.99e30_dp) then
			r(n) = 0.0_dp
			a(n) = 0.0_dp
		else
			r(n) = (-3.0_dp/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
			a(n) = 0.5
		end if
	
		call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
		
	end subroutine spline_real

	!
	! Returns a primitive function to a function y(t)exp(+i*w*t),
	! i.e. it calculates integral
    !
	!         t
	!        /
	!       |
	!       |dt' y(t')*exp(i*w*t')
	!       |
	!       /
	!      0
	!
	! for all values of t
	!
	subroutine primitive(x,y,w,py,D0)
		real(dp), dimension(:), intent(in) :: x,y
		real(dp), dimension(:), intent(in) :: w
		complex(dpc), optional, intent(in)     :: D0
		
		real(dp),     dimension(size(x)) :: d
		complex(dpc), dimension(size(x),size(w)), intent(out) :: py
		
		integer(i4b) :: m,i,k,n
		real(dp)      :: x1,x2,y1,y2,dd1,dd2
		real(dp), dimension(4) :: A
      	real(dp)      :: BeB,AeA,B2eB,A2eA,dI0,dI1,dI2,dI3
      	real(dp)      :: I0,I1,I2,I3 
       	complex(dpc) :: BeB_c,AeA_c,B2eB_c,A2eA_c,dI0_c,dI1_c,dI2_c,dI3_c
      	complex(dpc) :: I0_c,I1_c,I2_c,I3_c 
      	complex(dpc) :: eA,eB,iom
      	complex(dpc) :: der
		
		der = 0.0_dp
		if (present(D0)) der = D0
		
  		call spline_real(x,y,der,(0.0_dp,0.0_dp),d)
  
  		m = size(x) - 1
		n = size(w)
		
    	py = 0.0_dp
		
		do k = 1, n
  		if (abs(w(k)) < 1.0e-12_dp) then
  
    		!Ilast = 0
    		py(1,k) = 0.0_dp
    		do i = 1,m
  
      			x1 = x(i)
      			x2 = x(i+1)
      			y1 = y(i)
      			y2 = y(i+1)
      			dd1 = d(i)
     			dd2 = d(i+1)
  
      			call get_abcd(x1,x2,y1,y2,dd1,dd2,A)
  
      			BeB = x2*x2
      			AeA = x1*x1
      			B2eB = x2*BeB
      			A2eA = x1*AeA
      			dI0 = x2 - x1
      			dI1 = BeB - AeA
      			dI2 = B2eB - A2eA
      			dI3 = x2*B2eB - x1*A2eA
    
      			I0 = dI0
      			I1 = dI1/2.0_dp
      			I2 = dI2/3.0_dp
      			I3 = dI3/4.0_dp
  
      			py(i+1,k) = py(i,k) + A(1)*I0 + A(2)*I1 + A(3)*I2 + A(4)*I3;
  
    		end do
    		
  		else 
    
    		!print *, w(k)
    		py(1,k) = 0;
    		do i = 1,m
  
      			x1 = x(i)
      			x2 = x(i+1)
      			y1 = y(i)
      			y2 = y(i+1)
      			dd1 = d(i)
      			dd2 = d(i+1)
  
     	 		call get_abcd(x1,x2,y1,y2,dd1,dd2,A)
  
      			eA = exp((0.0_dp,1.0_dp)*w(k)*x1)
      			eB = exp((0.0_dp,1.0_dp)*w(k)*x2)
      			iom = (0.0_dp,1.0_dp)/w(k)
    
      			BeB_c = x2*eB
      			AeA_c = x1*eA
      			B2eB_c = x2*BeB_c
      			A2eA_c = x1*AeA_c
      			dI0_c = eB - eA
      			dI1_c = BeB_c - AeA_c
      			dI2_c = B2eB_c - A2eA_c
      			dI3_c = x2*B2eB_c - x1*A2eA_c
    
      			I0_c = -iom*dI0_c
      			I1_c = -iom*(dI1_c - I0_c)
      			I2_c = -iom*(dI2_c - 2.0_dp*I1_c)
      			I3_c = -iom*(dI3_c - 3.0_dp*I2_c)
    
      			py(i+1,k) = py(i,k) + A(1)*I0_c + A(2)*I1_c + &
      					A(3)*I2_c + A(4)*I3_c 
  
    		end do
  
  		end if
  		
  		end do

	end subroutine primitive

	subroutine primitive_cmplx(x,y,w,py,D0)
		real(dp), dimension(:), intent(in) :: x
		complex(dpc), dimension(:), intent(in) :: y
		real(dp), dimension(:), intent(in) :: w
		complex(dpc), optional, intent(in)     :: D0
		complex(dpc), dimension(size(x),size(w)), intent(out) :: py

		real(dp), dimension(size(x)) :: y_re, y_im
		complex(dpc), dimension(size(x),size(w)) :: py_re, py_im

		y_re = real(y)
		y_im = aimag(y)

		call primitive(x,y_re,w,py_re)
		call primitive(x,y_im,w,py_im)

		py = py_re + py_im * CMPLX(0,1)


	end subroutine primitive_cmplx
	
	
	!
	! Definite integral
	!
	function intsplin(x,y) result(res)
		real(dp), dimension(:), intenT(in) :: x,y
		!local
		real(dp), dimension(1) :: w
		complex(dpc), dimension(size(x),1) :: py
	    real(dp) :: res
	    
	    w(1) = 0.0_dp
	    call primitive(x,y,w,py)
		
		res = real(py(size(x),1))
		
	end function
	
	!
	! Cubic spline interpolation - complex
	!
	subroutine spline_cmplx(x,y,yp1,ypn,y2)
		real(dp), dimension(:), intent(in)  :: x
		complex(dpc), dimension(:), intent(in)  :: y
		complex(dpc), intent(in)                :: yp1, ypn
		complex(dp), dimension(:), intent(out) :: y2

		real(dp), dimension(:), allocatable :: re_y, im_y, re_y2, im_y2

		allocate(re_y2( size(y)) )
		allocate(im_y2( size(y)) )

		call spline_real(x,real(y),yp1,ypn,re_y2)
		call spline_real(x,aimag(y),yp1,ypn,im_y2)

		y2 = re_y2 + cmplx(0,1)*im_y2

		deallocate(re_y2)
		deallocate(im_y2)


	end subroutine spline_cmplx

	!
	! Cubic spline interpolation - complex
	!
	function splint_cmplx(xa,ya,y2a,x) result(ret)
		real(dp), dimension(:), intent(in) :: xa
		complex(dpc), dimension(:), intent(in) :: ya,y2a
		real(dp), intent(in) :: x
		complex(dpc) :: ret

		real(dp) :: A, B

		A = splint_real(xa,real(ya),real(y2a),x)
		B = splint_real(xa,aimag(ya),aimag(y2a),x)

		ret = A + cmplx(0,1)*B

	end function splint_cmplx

	!
	! Cubic spline interpolation for matrix
	!
	subroutine spline_matrix(x,y,yp1,ypn,y2)
		complex(dpc), dimension(:,:,:), intent(in)  	:: y
		real(dp), dimension(:), intent(in)  			:: x
		complex(dpc), intent(in), dimension(:,:)		:: yp1, ypn
		complex(dpc), dimension(:,:,:), intent(out) :: y2

		integer(i4b) :: i,j

		if(size(y,1) /= size(yp1,1) .or. size(y,2) /= size(yp1,2) .or. &
			size(y,1) /= size(ypn,1) .or. size(y,2) /= size(ypn,2) .or. size(x) /= size(y,3) .or.&
			size(y,1) /= size(y2,1) .or. size(y,2) /= size(y2,2) .or. size(y,3) /= size(y2,3) ) then

			call print_error_message(-1,"dimension error in splint_matrix")
		end if

		do i=1,size(y,1)
		do j=1,size(y,2)

			call spline(x,y(i,j,:),yp1(i,j),ypn(i,j),y2(i,j,:))

		end do
		end do

	end subroutine spline_matrix

	!
	! Cubic spline interpolation for matrix
	!
	subroutine splint_matrix(xa,ya,y2a,x,y)
		complex(dpc), dimension(:,:,:), intent(in) :: ya,y2a
		real(dp), dimension(:), intent(in) :: xa
		real(dp), intent(in) :: x
		complex(dpc), dimension(:,:), intent(out) :: y

		integer(i4b) :: i,j

		if(size(ya,1) /= size(y2a,1) .or. size(ya,2) /= size(y2a,2) .or. size(ya,3) /= size(y2a,3) .or. &
			size(xa) /= size(ya,3) .or. size(ya,1) /= size(y,1) .or. size(ya,2) /= size(y,2) ) then

			call print_error_message(-1,"dimension error in splint_matrix")
		end if

		do i=1,size(ya,1)
		do j=1,size(ya,2)

			y(i,j) = splint_cmplx(xa,ya(i,j,:),y2a(i,j,:),x)

		end do
		end do

	end subroutine splint_matrix

	
	
!***********************
!  Auxiliary functions
!***********************

	subroutine get_abcd(x1,x2,y1,y2,dd1,dd2,abcd_vec)
		real(dp), intent(in)                   :: x1,x2,y1,y2,dd1,dd2
		real(dp), dimension(4), intent(out)    :: abcd_vec	
		
		real(dp), dimension(4,4) :: M,M1
		real(dp), dimension(4)   :: vec
		
  		M(4,4)   = 0.0_dp
  		M(1:2,1) = 1.0_dp
  		M(3:4,2) = 1.0_dp
  		M(1,2) = x1
  		M(1,3) = M(1,2)*x1
  		M(1,4) = M(1,3)*x1
  		M(2,2) = x2
  		M(2,3) = M(2,2)*x2
  		M(2,4) = M(2,3)*x2
  		M(3,3) = 2.0_dp*x1
  		M(3,4) = 3.0_dp*x1*x1
  		M(4,3) = 2.0_dp*x2
  		M(4,4) = 3.0_dp*x2*x2
  
  		call inv(M,M1)
  
  		vec(1) = y1 ![y1,y2,dd1,dd2] !.'
  		vec(2) = y2
  		vec(3) = dd1
  		vec(4) = dd2

  		abcd_vec = matmul(M1,vec)

	end subroutine get_abcd


	subroutine tridag(a,b,c,r,u)
		real(dp), dimension(:), intent(in)  :: a,b,c,r
		real(dp), dimension(:), intent(out) :: u
		real(dp), dimension(size(b)) :: gam
		integer(i4b) :: n,j
		real(dp) :: bet
		
		n = size(b)
		bet = b(1)
		if (bet == 0.0_dp) stop "Error in tridag (1)"
		u(1) = r(1)/bet
		do j = 2,n
			gam(j) = c(j-1)/bet
			bet=b(j)-a(j-1)*gam(j)
			!print *, j,bet
			if (bet == 0.0_dp) stop "Error in tridag (2)"
			u(j) = (r(j)-a(j-1)*u(j-1))/bet
		end do
		do j = n-1,1,-1
			u(j) = u(j)-gam(j+1)*u(j+1)
		end do
		
	end subroutine tridag

	function locate(xx,x) result(loc)
		real(dp), dimension(:), intent(in) :: xx
		real(dp), intent(in) :: x
		integer(i4b) :: loc
		
		integer(i4b) :: n,j1,jm,ju
		logical :: ascnd
		
		n = size(xx)
		ascnd = (xx(n) >= xx(1))
		j1 = 0
		ju = n + 1
		do
			if (ju-j1 <= 1) exit
			jm = (ju+j1)/2
			if (ascnd .eqv. (x >= xx(jm))) then
				j1=jm
			else
				ju=jm
			end if
		end do
		if (x == xx(1)) then
			loc = 1
		else if (x == xx(n)) then
			loc = n - 1
		else
			loc = j1
		end if
					
	end function locate
		


end module numer_interp
