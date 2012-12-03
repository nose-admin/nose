!
! Random Number Generators
!
! - Uniform distribution
! - Normal distribution
! - Double-well potential distribution
!
!
! Tomas Mancal, 2006-11-23
!
! TODO: special random generators should go into some modules
!
module numer_random

	use std_types
    use std_io
	use numer_interp

	implicit none

	!
	! Public entities
	!

	! initiators
	public :: init_random
	public :: init_random_PotentialSingle
	public :: init_random_PotentialDouble

	! gets and sets
	public :: random_get_seedSize
	public :: random_get_seed
	public :: random_set_seed

	! Generators
	public :: random_Uniform
	public :: random_Normal
	public :: random_PotentialSingle
	public :: random_PotentialDouble

	! cleaners
	public :: clean_random

    public :: random_make_steps

	!
	! Private entities
	!
	private :: nseed, seeds
	private  :: mk1, mk2, mQ2, mE2, mJ, maa, mbeta, mZ
	private  :: xm, xw, MM
	private  :: seed_size, get_seed, set_seed
	private  :: Vtot, V1, V2, Venv, resety2, checky2
	private  :: clean, clean_into_file, init_from_file, make_steps

	integer :: nseed
	integer, dimension(:), allocatable :: seeds
	real(dp) ::	mk1, mk2, mQ2, mE2, mJ, maa, mbeta, mZ, msig
	real(dp) :: xm, xw, MM

	!
	! Interfaces
	!
	interface random_get_seedSize
		module procedure seed_size
	end interface random_get_seedSize
	interface random_get_seed
		module procedure get_seed
	end interface random_get_seed
	interface random_set_seed
		module procedure set_seed
	end interface random_set_seed
	interface init_random
		module procedure init, init_from_file, init_scratch
	end interface init_random
	interface clean_random
		module procedure clean, clean_into_file
	end interface clean_random
    interface random_make_steps
        module procedure make_steps
    end interface


contains

	function seed_size() result(res)
		integer(i4b):: res
		integer :: n
		call random_seed(n)
		res = n
	end function seed_size

	subroutine get_seed(seed)
		integer, dimension(:), intent(out) :: seed
		call random_seed(get=seed)
	end subroutine get_seed

	subroutine set_seed(seed)
		integer, dimension(:), intent(in) :: seed
		call random_seed(put=seed)
	end subroutine set_seed

    !
    ! Steps by n seeds forward
    !
    subroutine make_steps(n)
        integer, intent(in) :: n
        real, dimension(n) :: rr
        call random_number(rr)
    end subroutine make_steps


	!
	! Initialization of the Random Generator from an existing file
	! with the seed
	!
	! TODO: Check if the file exist, return errors
	! TODO: Check if unit 10 is free, select unit automatically
	!
	subroutine init_from_file(filename)
		character(len=*), intent(in) :: filename
		! local
		integer(i4b) :: i

      	open(10,file=filename)

      	read(10,*) nseed
      	allocate(seeds(nseed))
      	do i = 1, nseed
      		read(10,*) seeds(i)
      	end do

      	! initiation
      	call init(seeds)

      	close(10)

    end subroutine init_from_file

	!
	!
	!
	subroutine init(seed)
		integer, intent(in), dimension(:) :: seed
		real(dp) :: bkk1_T0, bkk2_T0, T0, x0, E2, J, T, kb, aa
		integer(i4b) :: n
		call print_log_message("Initialization of random module", 8)
		n = seed_size()
		if (n == size(seed)) then
			call set_seed(seed)
		else
			call print_error_message(-1,"init_random - wrong size of the seed array")
		end if

! 		bkk1_T0  =  150
! 		bkk2_T0  =  150
! 		T0       =  77
! 		x0      =   500
! 		E2      = 50
! 		J       = 80
! 		aa      = 1.0
! 		T        =  77
!
! 		print *, "calling init"
! 		call Init_RandPotentialDouble(bkk1_T0,bkk2_T0,T0,x0,E2,J,aa,T)
! 		print *, '..done'

		call print_log_message("... random module done.", 8)

	end subroutine init

    !
    ! Initializing with default seed
    !
    subroutine init_scratch()
        integer, dimension(:), allocatable :: seed
        integer :: n
        call random_seed(size=n)
        allocate(seed(n))
        call random_seed(get=seed)
        call init(seed)
    end subroutine init_scratch

	!
	! Initialization of the Double-Well Potential Disorder Generator
	!
	subroutine init_random_PotentialDouble(bkk1_T0,bkk2_T0,T0,x0,E2,J,aa,T)
		real(dp), intent(in) :: bkk1_T0, bkk2_T0, T0, x0, E2, J, aa, T

		! local
		real(dp) :: kb, bt
		real(dp) :: bt0, bkk1, bkk2, k1, k2, N
		real(dp), dimension(:), allocatable :: x
		real(dp) :: xmax, xmin, dx
		integer(i4b) :: nx, i, ccount, jj, Nk
		real(dp), dimension(:), allocatable :: y1,y2
		real(dp), dimension(:), allocatable :: zr
		real(dp) :: ConvET

        ConvET = 2.0_dp*acos(-1.0_dp)*2.99792458d-5

		!bkk1_T0  =  150
		!bkk2_T0  =  150
		!T0       =  77
		mQ2      = x0 ! 500
		mE2      = E2 ! 50
		mJ       = J  !80
		maa      = aa !1.0

		!T        =  77

		kb    =  0.6950387*ConvET

		mZ    = 1.0_dp
		mbeta = 1.0/(T*kb)
		bt = mbeta
		bt0   = 1.0/(T0*kb)
		bkk1  = bkk1_T0*sqrt(bt0/bt)
		bkk2  = bkk1_T0*sqrt(bt0/bt)
		mk1   = 1.0/(bt*(bkk1**2))
		mk2   = 1.0/(bt*(bkk2**2))
		N     = 1.0 + exp(-bt*E2)
		xm    = x0*(exp(-bt*E2)/N)
		xw    = (x0 + (1.0/N)*bkk1 + (exp(-bt*E2)/N)*bkk2)/2.0

		xmax = xm+4.0*xw
		xmin = xm-4.0*xw
		dx   = (xmax - xmin)/1000_dp
		nx   = int((xmax - xmin)/dx)

		allocate(x(nx),y1(nx),y2(nx))
		do i = 1, nx
			x(i) = xmin + (i-1)*dx
		end do
		y1   = exp(-bt*Vtot(x));

		mz   = intsplin(x,y1);
		y1   = y1/mz;

		! envelope distribution
		y2 = exp(-Venv(x,xm,xw))/(xw*sqrt(2.0*PI_D));

		! set M
		MM     = 1.0
		ccount = 0
		do
  			MM    = resety2(y1,y2,MM)
  			ccount = ccount + 1
  			if (checky2(y1,y2,MM) == 1) then
    			exit
  			end if
  			if (Ccount > 10) then
    			print *, "Didn't manage to converge after trying 10 times";
  				stop
  			end if
		end do

		print *, "Going to test"
		!
		! test the method
		!
! 		Nk = 100000
! 		allocate(zr(Nk))
!
! 		call RandPotentialDouble(zr)
!
! 		open(unit=10,file='dist.dat')
! 		do i = 1, Nk
! 			write(10,*) zr(i)
! 		end do
! 		close(10)

	end subroutine init_random_PotentialDouble

	!
	! Random thermal distribution on quadratic single-well potential
	!
	subroutine init_random_PotentialSingle(bkk_T0,T0,aa,T)
		real(dp), intent(in) :: bkk_T0, T0, aa, T
		! local
		real(dp) :: beta, bet0
		real(dp) :: x,u,fx,gx,M,r
		real(dp), dimension(1) :: xa
		real(dp) :: kb
		real(dp) :: ConvET

      	ConvET = 2.0_dp*acos(-1.0_dp)*2.99792458d-5
		kb = 0.6950387*ConvET

		bet0 = 1.0_dp/(kb*T0)
		beta = 1.0_dp/(kb*T)

		msig = bkk_T0*sqrt(bet0/beta)

	end subroutine init_random_PotentialSingle

	!
	! Uniform distributions
	!
	subroutine random_Uniform(res)
		real(dp), dimension(:), intent(out) :: res
		call random_number(res)
	end subroutine random_Uniform

	!
	! Normal Random Distribution
	!
	subroutine random_Normal(res)
		real(dp), dimension(:), intent(out) :: res
		real(dp) :: r1, r2,q1, qm, qp,pq,ee
		real(dp), dimension(2) :: rr

		integer(i4b) :: NN, i

		NN = size(res)

		i = 0
		do

			call random_number(rr)

			r1 = 2.0_dp*rr(1) - 1.0_dp
			r2 = 2.0_dp*rr(2) - 1.0_dp

			q1 = r1**2 + r2**2

			if ((q1 > 0.0_dp).and.(q1 <= 1.0_dp)) then
				i = i + 1
				res(i) = r1*sqrt(-2.0_dp*log(q1)/q1)
				if (i >= NN) exit
				i = i + 1
				res(i) = r2*sqrt(-2.0_dp*log(q1)/q1)
				if (i >= NN) exit
			end if

		end do

	end subroutine random_Normal

	!
	!
	!
	subroutine random_PotentialSingle(res)
		real(dp), dimension(:), intent(out) :: res

		call random_Normal(res)
		res  = msig*res

	end subroutine random_PotentialSingle


	!
	! Random Distribution on the double-well potential
	!
	subroutine random_PotentialDouble(res)
		real(dp), dimension(:), intent(out) :: res
		integer(i4b)                        :: jj, i
		real(dp), dimension(size(res))      :: xr, yr
		real(dp), dimension(1)              :: xx,yx,yz
		integer(i4b) :: Nk

		Nk = size(res);
		res(1:Nk) = 0.0
		jj = 0

		main: do

			call random_Normal(xr)
			call random_number(yr)
			xr = xw*xr + xm

			i = 0
			do

				i = i + 1

  				xx(1) = xr(i)
  				yx = exp(-mbeta*Vtot(xx))/mz
  				yz = exp(-Venv(xx,xm,xw))/(xw*sqrt(2.0*PI_D))
  				if ( yr(i) < yx(1)/(MM*yz(1)) ) then
    				jj = jj + 1
    				res(jj) = xx(1)
    				if (jj >= Nk) exit MAIN
  				end if

  				if (i == Nk) exit

			end do

		end do main

	end subroutine random_PotentialDouble


	!
	! Saving seed and cleaning
	!
	subroutine clean_into_file(filename)
		character(len=*), intent(in) :: filename
		! local
		integer(i4b) :: i

      	open(10,file=filename)

      	write(10,*) nseed
      	do i = 1, nseed
      		write(10,*) seeds(i)
      	end do

    	call clean

      	close(10)

	end subroutine clean_into_file

	!
	! Cleaning the module
	!
	subroutine clean

	    deallocate(seeds)

	end subroutine clean


!
! Auxiliarly functions
!

	!
	! Well tot
	!
	function Vtot(q) result (res)
		real(dp), intent(in), dimension(:) :: q
		real(dp), dimension(size(q))       :: res

		res = (V1(q)+V2(q) - sqrt((V1(q)-V2(q))**2 + 4.0_dp*(mJ**2)))/2.0_dp

	end function Vtot


	!
	! Well 1
	!
	function V1(q) result (res)
		real(dp), intent(in), dimension(:) :: q
		real(dp), dimension(size(q))       :: res

		res = mk1*(q**2)/2.0_dp

	end function V1

	!
	! Well 2
	!
	function V2(q) result (res)
		real(dp), intent(in), dimension(:) :: q
		real(dp), dimension(size(q))       :: res

		res = mk2*((q-mQ2)**2)/2.0_dp + mE2

	end function V2

	!
	! Envelop distribution (Normal Gaussian)
	!
	function Venv(x,xm,xw) result (res)
		real(dp), intent(in) :: xm,xw
		real(dp), dimension(:) :: x
		real(dp), dimension(size(x)) :: res

  		res = ((x-xm)**2)/(2*(xw**2))

	end function Venv

	!
	! Calculates the multiplication constant M
	!
	function resety2(y1, y2, M_in) result (M)
		real(dp), intent(in)               :: M_in
		real(dp), intent(in), dimension(:) :: y1, y2
		real(dp) :: M
		! local
		real(dp), dimension(size(y1))      :: yy
		integer(i4b), dimension(1) :: ind
		real(dp)     :: val, val1, val2, M_new

  		yy     = (M_in*y2 - y1)
  		val    = maxval(-yy)
  		ind    = maxloc(-yy)
  		val1   = y1(ind(1))
  		val2   = M_in*y2(ind(1))
  		M_new  = val1/val2
  		M      = M_new*M_in

	end function resety2

	!
	! Check if the envelop works
	!
	function checky2(y1,y2,M) result (ss)
		real(dp), intent(in)               :: M
		real(dp), intent(in), dimension(:) :: y1, y2
		! local
		integer(i4b) :: ss
		integer(i4b), dimension(1) ::indm
		real(dp)     :: ym
		real(dp), dimension(size(y1)) :: yy

  		yy   = M*y2 - y1;
  		ss   = 1;
  		ym   = minval(yy);
  		indm = minloc(yy);
  		if (ym < 0.0) then
    		ss = -1;
  		end if

	end function


















!
!  I
!  I   To be delete
! \ /
!
!********** BEGIN of Random Gauss Function Unit in Mino's program
!
!   IMSL ROUTINE NAME   - GGNQF
!
!-----------------------------------------------------------------------
!
!   LATEST REVISION     - NOVEMBER 1, 1984
!
!   PURPOSE             - NORMAL RANDOM DEVIATE GENERATOR - FUNCTION
!                           FORM OF GGNML
!
!   USAGE               - FUNCTION GGNQF (DSEED)
!
!   ARGUMENTS    GGNQF  - RESULTANT NORMAL (0,1) DEVIATE.
!                DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
!                           ASSIGNED AN INTEGER VALUE IN THE
!                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
!                           DSEED IS REPLACED BY A NEW VALUE TO BE
!                           USED IN A SUBSEQUENT CALL.
!
!   REQD. IMSL ROUTINES - MDNRIS,MERFI,UERTST,UGETIO
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                          CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!-----------------------------------------------------------------------
!
      REAL FUNCTION RandGauss (DSEED)
!                                  SPECIFICATIONS FOR ARGUMENTS
      DOUBLE PRECISION   DSEED
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IER
      REAL               GGNQFX,XGGNQF
      DOUBLE PRECISION   D2P31M,D2PN31
!                                  D2P31M = (2**31) - 1
!                                  D2PN31 = (2**31)
      DATA               D2P31M/2147483647.D0/
      DATA               D2PN31/2147483648.D0/
!                                  GENERATE A RANDOM (0,1) DEVIATE.
!                                  16807 = (7**5)
!                                  FIRST EXECUTABLE STATEMENT
      DSEED = DMOD(16807.D0*DSEED,D2P31M)
      GGNQFX = DSEED / D2PN31
!                                  TRANSFORM TO NORMAL DEVIATE
      CALL MDNRIS(GGNQFX,XGGNQF,IER)
      RandGauss = XGGNQF

      END FUNCTION

!   IMSL ROUTINE NAME   - MDNRIS
!
!-----------------------------------------------------------------------
!
!   LATEST REVISION     - JUNE 1, 1981
!
!   PURPOSE             - INVERSE STANDARD NORMAL (GAUSSIAN)
!                           PROBABILITY DISTRIBUTION FUNCTION
!
!   USAGE               - CALL MDNRIS (P,Y,IER)
!
!   ARGUMENTS    P      - INPUT VALUE IN THE EXCLUSIVE RANGE (0.0,1.0)
!                Y      - OUTPUT VALUE OF THE INVERSE NORMAL (0,1)
!                           PROBABILITY DISTRIBUTION FUNCTION
!                IER    - ERROR PARAMETER (OUTPUT)
!                         TERMINAL ERROR
!                           IER = 129 INDICATES P LIES OUTSIDE THE LEGAL
!                             RANGE. PLUS OR MINUS MACHINE INFINITY IS
!                             GIVEN AS THE RESULT (SIGN IS THE SIGN OF
!                             THE FUNCTION VALUE OF THE NEAREST LEGAL
!                             ARGUMENT).
!
!   REQD. IMSL ROUTINES - MERFI,UERTST,UGETIO
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!-----------------------------------------------------------------------
!
      SUBROUTINE MDNRIS (P,Y,IER)
!                                  SPECIFICATIONS FOR ARGUMENTS
      REAL               P,Y
      INTEGER            IER
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               EPS,G0,G1,G2,G3,H0,H1,H2,A,W,WI,SN,SD
      REAL               SIGMA,SQRT2,X,XINF
!      DATA               XINF/Z7FFFFFFF/
      DATA		 XINF/1.E30/
      DATA               SQRT2/1.414214/
!      DATA               EPS/Z3C100000/
      DATA		 EPS/1.E-30/
      DATA               G0/.1851159E-3/,G1/-.2028152E-2/
      DATA               G2/-.1498384/,G3/.1078639E-1/
      DATA               H0/.9952975E-1/,H1/.5211733/
      DATA               H2/-.6888301E-1/
!                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (P .GT. 0.0 .AND. P .LT. 1.0) GO TO 5
      IER = 129
      SIGMA = SIGN(1.0,P)
      Y = SIGMA * XINF
      GO TO 9000
    5 IF(P.LE.EPS) GO TO 10
      X = 1.0 -(P + P)
      CALL MERFI (X,Y,IER)
      Y = -SQRT2 * Y
      GO TO 9005
!                                  P TOO SMALL, COMPUTE Y DIRECTLY
   10 A = P+P
      W = SQRT(-ALOG(A+(A-A*A)))
!                                  USE A RATIONAL FUNCTION IN 1./W
      WI = 1./W
      SN = ((G3*WI+G2)*WI+G1)*WI
      SD = ((WI+H2)*WI+H1)*WI+H0
      Y = W + W*(G0+SN/SD)
      Y = -Y*SQRT2
      GO TO 9005
 9000 CONTINUE
!      CALL UERTST(IER,6HMDNRIS)
 9005 RETURN
      END SUBROUTINE

! C   IMSL ROUTINE NAME   - MERFI
! C
! C-----------------------------------------------------------------------
! C
! C   LATEST REVISION     - JANUARY 1, 1978
! C
! C   PURPOSE             - INVERSE ERROR FUNCTION
! C
! C   USAGE               - CALL MERFI (P,Y,IER)
! C
! C   ARGUMENTS    P      - INPUT VALUE IN THE EXCLUSIVE RANGE (-1.0,1.0)
! C                Y      - OUTPUT VALUE OF THE INVERSE ERROR FUNCTION
! C                IER    - ERROR PARAMETER (OUTPUT)
! C                         TERMINAL ERROR
! C                           IER = 129 INDICATES P LIES OUTSIDE THE LEGAL
! C                             RANGE. PLUS OR MINUS MACHINE INFINITY IS
! C                             GIVEN AS THE RESULT (SIGN IS THE SIGN OF
! C                             THE FUNCTION VALUE OF THE NEAREST LEGAL
! C                             ARGUMENT).
! C
! C   REQD. IMSL ROUTINES - UERTST,UGETIO
! C
! C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
! C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
! C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
! C
! C-----------------------------------------------------------------------
! C
      SUBROUTINE MERFI (P,Y,IER)
!C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL               P,Y
      INTEGER            IER
!C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               A,B,X,Z,W,WI,SN,SD,F,Z2,RINFM,A1,A2,A3,B0,B1, &
                        B2,B3,C0,C1,C2,C3,D0,D1,D2,E0,E1,E2,E3,F0,F1, &
                        F2,G0,G1,G2,G3,H0,H1,H2,SIGMA
      DATA               A1/-.5751703/,A2/-1.896513/,A3/-.5496261E-1/
      DATA               B0/-.1137730/,B1/-3.293474/,B2/-2.374996/
      DATA               B3/-1.187515/
      DATA               C0/-.1146666/,C1/-.1314774/,C2/-.2368201/
      DATA               C3/.5073975E-1/
      DATA               D0/-44.27977/,D1/21.98546/,D2/-7.586103/
      DATA               E0/-.5668422E-1/,E1/.3937021/,E2/-.3166501/
      DATA               E3/.6208963E-1/
      DATA               F0/-6.266786/,F1/4.666263/,F2/-2.962883/
      DATA               G0/.1851159E-3/,G1/-.2028152E-2/
      DATA               G2/-.1498384/,G3/.1078639E-1/
      DATA               H0/.9952975E-1/,H1/.5211733/
      DATA               H2/-.6888301E-1/
!      DATA               RINFM/Z7FFFFFFF/
      DATA               RINFM/1.E30/
!C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      X = P
      SIGMA = SIGN(1.0,X)
!C                                  TEST FOR INVALID ARGUMENT
      IF (.NOT.(X.GT.-1. .AND. X.LT.1.)) GO TO 30
      Z = ABS(X)
      IF (Z.LE. .85) GO TO 20
      A = 1.-Z
      B = Z
!C                                  REDUCED ARGUMENT IS IN (.85,1.),
!C                                     OBTAIN THE TRANSFORMED VARIABLE
    5 W = SQRT(-ALOG(A+A*B))
      IF (W.LT.2.5) GO TO 15
      IF (W.LT.4.) GO TO 10
!C                                  W GREATER THAN 4., APPROX. F BY A
!C                                     RATIONAL FUNCTION IN 1./W
      WI = 1./W
      SN = ((G3*WI+G2)*WI+G1)*WI
      SD = ((WI+H2)*WI+H1)*WI+H0
      F = W + W*(G0+SN/SD)
      GO TO 25
!C                                  W BETWEEN 2.5 AND 4., APPROX. F
!C                                     BY A RATIONAL FUNCTION IN W
   10 SN = ((E3*W+E2)*W+E1)*W
      SD = ((W+F2)*W+F1)*W+F0
      F = W + W*(E0+SN/SD)
      GO TO 25
!C                                  W BETWEEN 1.13222 AND 2.5, APPROX.
!C                                     F BY A RATIONAL FUNCTION IN W
   15 SN = ((C3*W+C2)*W+C1)*W
      SD = ((W+D2)*W+D1)*W+D0
      F = W + W*(C0+SN/SD)
      GO TO 25
!C                                  Z BETWEEN 0. AND .85, APPROX. F
!C                                     BY A RATIONAL FUNCTION IN Z
   20 Z2 = Z*Z
      F = Z+Z*(B0+A1*Z2/(B1+Z2+A2/(B2+Z2+A3/(B3+Z2))))
!C                                  FORM THE SOLUTION BY MULT. F BY
!C                                     THE PROPER SIGN
   25 Y = SIGMA*F
      IER = 0
      GO TO 9005
!C                                  ERROR EXIT. SET SOLUTION TO PLUS
!C                                     (OR MINUS) INFINITY
   30 IER = 129
      Y = SIGMA * RINFM
 9000 CONTINUE
!c      CALL UERTST(IER,6HMERFI )
 9005 RETURN
      END SUBROUTINE

!c*************** END OF Random Gauss Function Unit **********



end module numer_random
