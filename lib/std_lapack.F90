!
! Simple interface to lapack based on original LAPACK95 package
! See: LAPACK95 Users' Guide, available at http://www.netlib.org/lapack95/lug95/
!
! The routines std_lapack_geev and std_lapack_gesv are more or less cut and pasted versions of
! la_geev and la_gesv subrotines of that package. The LAPACK95 licence seems to allow this.
!
! std_types are replacing the LAPACK95 types
! subroutine LSAME and ERINFO are used in their original version
!
! Author: Tomas Mancal, mancal@karlov.mff.cuni.cz
! Date: 2008/02/15
!   
module std_lapack
 
    use std_types
    
    IMPLICIT NONE

    ! Interface to lapack routines
    
    INTERFACE LA_GEEV_m

       SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, &
                        LDVR, WORK, LWORK, INFO )
         use std_types
         CHARACTER(LEN=1), INTENT(IN) :: JOBVL, JOBVR
         INTEGER, INTENT(IN) :: N, LDA, LDVL, LDVR, LWORK
         INTEGER, INTENT(OUT) :: INFO
         REAL(DP), INTENT(INOUT) :: A(LDA,*)
         REAL(DP), INTENT(OUT) :: VL(LDVL,*), VR(LDVR,*), WR(*), WI(*), &
                                 WORK(*)
      END SUBROUTINE DGEEV

       SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR,&
                        WORK, LWORK, RWORK, INFO )
                        use std_types
         CHARACTER(LEN=1), INTENT(IN) :: JOBVL, JOBVR
         INTEGER, INTENT(IN) :: N, LDA, LDVL, LDVR, LWORK
         INTEGER, INTENT(OUT) :: INFO
         REAL(DPC), INTENT(OUT) :: RWORK(*)
         COMPLEX(DPC), INTENT(INOUT) :: A(LDA,*)
         COMPLEX(DPC), INTENT(OUT) :: VL(LDVL,*), VR(LDVR,*), W(*),      &
                                    WORK(*)
      END SUBROUTINE ZGEEV

    END INTERFACE

    INTERFACE LA_GESV_m

       SUBROUTINE DGESV( N, NRHS, A, LDA, PIV, B, LDB, INFO )
                        use std_types
         INTEGER, INTENT(IN) :: LDA, LDB, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: PIV(*)
         REAL(DP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
      END SUBROUTINE DGESV

       SUBROUTINE ZGESV( N, NRHS, A, LDA, PIV, B, LDB, INFO )
                        use std_types
         INTEGER, INTENT(IN) :: LDA, LDB, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: PIV(*)
         COMPLEX(DPC), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
      END SUBROUTINE ZGESV

    END INTERFACE

    INTERFACE std_lapack_gesv
    	MODULE PROCEDURE std_lapack_gesv_cmplx
    	MODULE PROCEDURE std_lapack_gesv_real
    END INTERFACE

    INTERFACE std_lapack_geev
    	MODULE PROCEDURE std_lapack_geev_cmplx
    	MODULE PROCEDURE std_lapack_geev_real
    END INTERFACE

contains


!    subroutine std_lapack_geev(AC,W,W2,S1,info)
!        real(dp), dimension(:,:), intent(out) :: S1
!        real(dp), dimension(:) :: W,W2
!        real(dp), dimension(:,:):: AC
!        integer :: info

!        call la_geev_i(AC,W,W2,VR=S1,INFO=info)
!    end subroutine std_lapack_geev    
    
!    subroutine std_lapack_gesv(Ai,B,info)
!        real(dp), dimension(:,:), intent(out) :: Ai
!        real(dp), dimension(:,:) :: B
!        integer :: info
!        call la_gesv_i(Ai,B,INFO=info)
!    end subroutine std_lapack_gesv    


    subroutine std_lapack_geev_cmplx( A, W, W2, VL, VR, INFO )
       IMPLICIT NONE
    !  .. SCALAR ARGUMENTS ..
       INTEGER, INTENT(OUT), OPTIONAL :: INFO
    !  .. ARRAY ARGUMENTS ..
       COMPLEX(dpc), INTENT(INOUT) :: A(:,:)
       COMPLEX(dpc), INTENT(OUT) :: W(:), W2(:)
       COMPLEX(dpc), INTENT(OUT), OPTIONAL, TARGET :: VL(:,:), VR(:,:)

    !  .. LOCAL PARAMETERS ..
       CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_GEEV'
    !  .. LOCAL SCALARS ..
       CHARACTER(LEN=1) :: LJOBVL, LJOBVR
       INTEGER, SAVE :: LWORK = 0
       INTEGER :: N, NN, LINFO, LD, ISTAT, ISTAT1, S1VL, S2VL, S1VR, S2VR
    !  .. LOCAL ARRAYS ..
       COMPLEX(dpc), TARGET :: LLVL(1,1), LLVR(1,1)
       COMPLEX(dpc), POINTER :: WORK(:)
       REAL(dp), POINTER :: RWORK(:) !!!!!!!!!!!!!!!!
    !  .. INTRINSIC FUNCTIONS ..
       INTRINSIC MAX, PRESENT, SIZE
    !  .. EXECUTABLE STATEMENTS ..
       LINFO = 0; ISTAT = 0; N = SIZE(A,1); LD = MAX(1,N)
       IF( PRESENT(VL) )THEN; S1VL = SIZE(VL,1); S2VL = SIZE(VL,2); LJOBVL = 'V'
       ELSE; S1VL = 1; S2VL = 1; LJOBVL = 'N'; END IF
       IF( PRESENT(VR) )THEN; S1VR = SIZE(VR,1); S2VR = SIZE(VR,2); LJOBVR = 'V'
       ELSE; S1VR = 1; S2VR = 1; LJOBVR = 'N'; END IF
    !  .. TEST THE ARGUMENTS
       IF( N < 0 .OR. SIZE(A,2) /= N )THEN; LINFO = -1
       !HACK
       ELSE IF( SIZE( W ) /= N )THEN; LINFO = -2
       !ELSE IF( SIZE( WI ) /= N )THEN; LINFO = -3
       ELSE IF( PRESENT(VL) .AND. ( S1VL /= N .OR. S2VL /= N ) )THEN; LINFO = -4
       ELSE IF( PRESENT(VR) .AND. ( S1VR /= N .OR. S2VR /= N ) )THEN; LINFO = -5
       ELSE IF( N > 0 )THEN
          ALLOCATE(RWORK(2*N))
          NN = 3; IF( LSAME(LJOBVL,'V').OR.LSAME(LJOBVR,'V') ) NN = NN + 1
          LWORK = MAX( 1, NN*N, LWORK); ALLOCATE(WORK(LWORK), STAT=ISTAT)
             IF( ISTAT /= 0 )THEN; DEALLOCATE(WORK,STAT=ISTAT1)
                LWORK = MAX( 1, NN*N ); ALLOCATE(WORK(LWORK), STAT=ISTAT)
                IF( ISTAT == 0) CALL ERINFO( -200, SRNAME, LINFO )
             END IF
          IF( ISTAT == 0 ) THEN
            IF( PRESENT(VL) )THEN
               IF( PRESENT(VR) )THEN
                   CALL la_GEEV_m( LJOBVL, LJOBVR, N, A, LD, W, &
                            VL, S1VL, VR, S1VR, WORK, LWORK, RWORK, LINFO )
           ELSE
              CALL la_GEEV_m( LJOBVL, LJOBVR, N, A, LD, W, &
                            VL, S1VL, LLVR, S1VR, WORK, LWORK, RWORK, LINFO )
           ENDIF
         ELSE
           IF( PRESENT(VR) )THEN
               CALL la_GEEV_m( LJOBVL, LJOBVR, N, A, LD, W, &
                            LLVL, S1VL, VR, S1VR, WORK, LWORK, RWORK, LINFO )
           ELSE
               CALL la_GEEV_m( LJOBVL, LJOBVR, N, A, LD, W, &
                            LLVL, S1VL, LLVR, S1VR, WORK, LWORK, RWORK, LINFO )
           ENDIF
         ENDIF
             IF( LINFO == 0 ) LWORK = INT(WORK(1)+1)
          ELSE; LINFO = -100; ENDIF
          DEALLOCATE(WORK, STAT=ISTAT1)
       ENDIF
       CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
    end subroutine std_lapack_geev_cmplx


    subroutine std_lapack_geev_real( A, WR, WI, VL, VR, INFO )
       IMPLICIT NONE
    !  .. SCALAR ARGUMENTS ..
       INTEGER, INTENT(OUT), OPTIONAL :: INFO
    !  .. ARRAY ARGUMENTS ..
       REAL(dp), INTENT(INOUT) :: A(:,:)
       REAL(dp), INTENT(OUT) :: WR(:), WI(:)
       REAL(dp), INTENT(OUT), OPTIONAL, TARGET :: VL(:,:), VR(:,:)
    
    !  .. LOCAL PARAMETERS ..
       CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_GEEV'
    !  .. LOCAL SCALARS ..
       CHARACTER(LEN=1) :: LJOBVL, LJOBVR
       INTEGER, SAVE :: LWORK = 0
       INTEGER :: N, NN, LINFO, LD, ISTAT, ISTAT1, S1VL, S2VL, S1VR, S2VR
    !  .. LOCAL ARRAYS ..
       REAL(dp), TARGET :: LLVL(1,1), LLVR(1,1)
       REAL(dp), POINTER :: WORK(:)
    !  .. INTRINSIC FUNCTIONS ..
       INTRINSIC MAX, PRESENT, SIZE
    !  .. EXECUTABLE STATEMENTS ..
       LINFO = 0; ISTAT = 0; N = SIZE(A,1); LD = MAX(1,N)
       IF( PRESENT(VL) )THEN; S1VL = SIZE(VL,1); S2VL = SIZE(VL,2); LJOBVL = 'V'
       ELSE; S1VL = 1; S2VL = 1; LJOBVL = 'N'; END IF
       IF( PRESENT(VR) )THEN; S1VR = SIZE(VR,1); S2VR = SIZE(VR,2); LJOBVR = 'V'
       ELSE; S1VR = 1; S2VR = 1; LJOBVR = 'N'; END IF
    !  .. TEST THE ARGUMENTS
       IF( N < 0 .OR. SIZE(A,2) /= N )THEN; LINFO = -1
       ELSE IF( SIZE( WR ) /= N )THEN; LINFO = -2
       ELSE IF( SIZE( WI ) /= N )THEN; LINFO = -3
       ELSE IF( PRESENT(VL) .AND. ( S1VL /= N .OR. S2VL /= N ) )THEN; LINFO = -4
       ELSE IF( PRESENT(VR) .AND. ( S1VR /= N .OR. S2VR /= N ) )THEN; LINFO = -5
       ELSE IF( N > 0 )THEN
          NN = 3; IF( LSAME(LJOBVL,'V').OR.LSAME(LJOBVR,'V') ) NN = NN + 1
          LWORK = MAX( 1, NN*N, LWORK); ALLOCATE(WORK(LWORK), STAT=ISTAT)
             IF( ISTAT /= 0 )THEN; DEALLOCATE(WORK,STAT=ISTAT1)
                LWORK = MAX( 1, NN*N ); ALLOCATE(WORK(LWORK), STAT=ISTAT)
                IF( ISTAT == 0) CALL ERINFO( -200, SRNAME, LINFO )
             END IF
          IF( ISTAT == 0 ) THEN
            IF( PRESENT(VL) )THEN
               IF( PRESENT(VR) )THEN
                   CALL la_GEEV_m( LJOBVL, LJOBVR, N, A, LD, WR, WI, &
                            VL, S1VL, VR, S1VR, WORK, LWORK, LINFO )
           ELSE
              CALL la_GEEV_m( LJOBVL, LJOBVR, N, A, LD, WR, WI, &
                            VL, S1VL, LLVR, S1VR, WORK, LWORK, LINFO )
           ENDIF
         ELSE
           IF( PRESENT(VR) )THEN
               CALL la_GEEV_m( LJOBVL, LJOBVR, N, A, LD, WR, WI, &
                            LLVL, S1VL, VR, S1VR, WORK, LWORK, LINFO )
           ELSE
               CALL la_GEEV_m( LJOBVL, LJOBVR, N, A, LD, WR, WI, &
                            LLVL, S1VL, LLVR, S1VR, WORK, LWORK, LINFO )
           ENDIF
         ENDIF  
             IF( LINFO == 0 ) LWORK = INT(WORK(1)+1)
          ELSE; LINFO = -100; ENDIF
          DEALLOCATE(WORK, STAT=ISTAT1)
       ENDIF
       CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
    end subroutine std_lapack_geev_real
    
    
    subroutine std_lapack_gesv_real( A, B, IPIV, INFO )
      IMPLICIT NONE
!     .. "Scalar Arguments" ..
      INTEGER, INTENT(OUT), OPTIONAL :: INFO
!     .. "Array Arguments" ..
      INTEGER, INTENT(OUT), OPTIONAL, TARGET :: IPIV(:)
      REAL(dp), INTENT(INOUT) :: A(:,:), B(:,:)

!     .. "Parameters" ..
      CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_GESV'
!     .. LOCAL SCALARS ..
      INTEGER :: LINFO, ISTAT, ISTAT1, SIPIV
!     .. "Local Pointers" ..
      INTEGER, POINTER :: LPIV(:)
!     .. "Intrinsic Functions" ..
      INTRINSIC SIZE, PRESENT, MAX
!     .. "Executable Statements" ..
      LINFO = 0; ISTAT = 0
      IF( PRESENT(IPIV) )THEN
         SIPIV = SIZE(IPIV)
      ELSE
         SIPIV = SIZE(A,1)
      END IF
!     .. "Test the arguments" ..
      IF( SIZE( A, 2 ) /= SIZE(A,1) .OR. SIZE(A,1) < 0 ) THEN
         LINFO = -1
      ELSE IF( SIZE( B, 1 ) /= SIZE(A,1) .OR. SIZE(B,2) < 0 ) THEN
         LINFO = -2
      ELSE IF( SIPIV /= SIZE(A,1) )THEN
            LINFO = -3
      ELSE IF ( SIZE(A,1) > 0 ) THEN
         IF( PRESENT(IPIV) )THEN
            LPIV => IPIV
         ELSE
            ALLOCATE( LPIV(SIZE(A,1)), STAT = ISTAT )
         END IF
         IF ( ISTAT == 0 ) THEN
!        .. "Call LAPACK77 routine" ..
            CALL la_GESV_m( SIZE(A,1), SIZE(B,2), A, MAX(1,SIZE(A,1)), LPIV, B, MAX(1,SIZE(A,1)), &
                           LINFO )
         ELSE
            LINFO = -100
         END IF
         IF( .NOT.PRESENT(IPIV) )THEN
            DEALLOCATE(LPIV, STAT = ISTAT1 )
         END IF
      END IF
      CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )

    end subroutine std_lapack_gesv_real

    subroutine std_lapack_gesv_cmplx( A, B, IPIV, INFO )
      IMPLICIT NONE
!     .. "Scalar Arguments" ..
      INTEGER, INTENT(OUT), OPTIONAL :: INFO
!     .. "Array Arguments" ..
      INTEGER, INTENT(OUT), OPTIONAL, TARGET :: IPIV(:)
      COMPLEX(dpc), INTENT(INOUT) :: A(:,:), B(:,:)

!     .. "Parameters" ..
      CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_GESV'
!     .. LOCAL SCALARS ..
      INTEGER :: LINFO, ISTAT, ISTAT1, SIPIV
!     .. "Local Pointers" ..
      INTEGER, POINTER :: LPIV(:)
!     .. "Intrinsic Functions" ..
      INTRINSIC SIZE, PRESENT, MAX
!     .. "Executable Statements" ..
      LINFO = 0; ISTAT = 0
      IF( PRESENT(IPIV) )THEN
         SIPIV = SIZE(IPIV)
      ELSE
         SIPIV = SIZE(A,1)
      END IF
!     .. "Test the arguments" ..
      IF( SIZE( A, 2 ) /= SIZE(A,1) .OR. SIZE(A,1) < 0 ) THEN
         LINFO = -1
      ELSE IF( SIZE( B, 1 ) /= SIZE(A,1) .OR. SIZE(B,2) < 0 ) THEN
         LINFO = -2
      ELSE IF( SIPIV /= SIZE(A,1) )THEN
            LINFO = -3
      ELSE IF ( SIZE(A,1) > 0 ) THEN
         IF( PRESENT(IPIV) )THEN
            LPIV => IPIV
         ELSE
            ALLOCATE( LPIV(SIZE(A,1)), STAT = ISTAT )
         END IF
         IF ( ISTAT == 0 ) THEN
!        .. "Call LAPACK77 routine" ..
            CALL la_GESV_m( SIZE(A,1), SIZE(B,2), A, MAX(1,SIZE(A,1)), LPIV, B, MAX(1,SIZE(A,1)), &
                           LINFO )
         ELSE
            LINFO = -100
         END IF
         IF( .NOT.PRESENT(IPIV) )THEN
            DEALLOCATE(LPIV, STAT = ISTAT1 )
         END IF
      END IF
      CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )

    end subroutine std_lapack_gesv_cmplx


     LOGICAL FUNCTION LSAME( CA, CB )
!
!  PURPOSE
!  =======
!
!  LSAME  TESTS IF CA IS THE SAME LETTER AS CB REGARDLESS OF CASE.
!
!  PARAMETERS
!  ==========
!
!  CA      (INPUT) CHARACTER*1
!  CB      (INPUT) CHARACTER*1
!          CHARACTERS TO BE COMPARED.
!
!  .. SCALAR ARGUMENTS ..
      CHARACTER*1, INTENT(IN) :: CA, CB
!  .. PARAMETERS ..
      INTEGER, PARAMETER      :: IOFF=32
!  .. LOCAL SCALARS ..
      INTEGER                 :: INTA, INTB, ZCODE
!  .. INTRINSIC FUNCTIONS ..
      INTRINSIC                  ICHAR
!
!  .. EXECUTABLE STATEMENTS ..
!
!  TEST IF THE CHARACTERS ARE EQUAL
!
      LSAME = CA == CB
!
!  NOW TEST FOR EQUIVALENCE
!
      IF( .NOT.LSAME )THEN
!
!     USE 'Z' RATHER THAN 'A' SO THAT ASCII CAN BE DETECTED ON PRIME
!     MACHINES, ON WHICH ICHAR RETURNS A VALUE WITH BIT 8 SET.
!     ICHAR('A') ON PRIME MACHINES RETURNS 193 WHICH IS THE SAME AS
!     ICHAR('A') ON AN EBCDIC MACHINE.
!
         ZCODE = ICHAR( 'Z' )
!
         INTA = ICHAR( CA )
         INTB = ICHAR( CB )
!
         IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 )THEN
!
!        ASCII IS ASSUMED - ZCODE IS THE ASCII CODE OF EITHER LOWER OR
!        UPPER CASE 'Z'.
!
            IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
            IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
!
         ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 )THEN
!
!        EBCDIC IS ASSUMED - ZCODE IS THE EBCDIC CODE OF EITHER LOWER OR
!        UPPER CASE 'Z'.
!
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.                         &
!    &       INTA.GE.145 .AND. INTA.LE.153 .OR.                         &
     &       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.                         &
     &       INTB.GE.145 .AND. INTB.LE.153 .OR.                         &
     &       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
!
         ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 )THEN
!
!        ASCII IS ASSUMED, ON PRIME MACHINES - ZCODE IS THE ASCII CODE
!        PLUS 128 OF EITHER LOWER OR UPPER CASE 'Z'.
!
            IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
         ENDIF
         LSAME = INTA == INTB
      ENDIF
      END FUNCTION LSAME

      SUBROUTINE ERINFO(LINFO, SRNAME, INFO, ISTAT)
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. IMPLICIT STATEMENT ..
         IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
         CHARACTER( LEN = * ), INTENT(IN)              :: SRNAME
         INTEGER             , INTENT(IN)              :: LINFO
         INTEGER             , INTENT(OUT), OPTIONAL   :: INFO
         INTEGER             , INTENT(IN), OPTIONAL    :: ISTAT
!  .. EXECUTABLE STATEMENTS ..
!         IF( ( LINFO < 0 .AND. LINFO > -200 ) .OR.                     &
!    &       ( LINFO > 0 .AND. .NOT.PRESENT(INFO) ) )THEN
      IF( ( ( LINFO < 0 .AND. LINFO > -200 ) .OR. LINFO > 0 )           &
     &           .AND. .NOT.PRESENT(INFO) )THEN
        WRITE (*,*) 'Program terminated in LAPACK95 subroutine ',SRNAME
        WRITE (*,*) 'Error indicator, INFO = ',LINFO
        IF( PRESENT(ISTAT) )THEN
          IF( ISTAT /= 0 ) THEN
            IF( LINFO == -100 )THEN
              WRITE (*,*) 'The statement ALLOCATE causes STATUS = ',    &
     &                    ISTAT
            ELSE
              WRITE (*,*) 'LINFO = ', LINFO, ' not expected'
            END IF
          END IF   
        END IF
        STOP
         ELSE IF( LINFO <= -200 ) THEN
           WRITE(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++'
           WRITE(*,*) '*** WARNING, INFO = ', LINFO, ' WARNING ***'
           IF( LINFO == -200 )THEN
             WRITE(*,*)                                                 &
     &        'Could not allocate sufficient workspace for the optimum'
             WRITE(*,*)                                                 &
     &        'blocksize, hence the routine may not have performed as'
             WRITE(*,*) 'efficiently as possible'
         ELSE
           WRITE(*,*) 'Unexpected warning'
         END IF
           WRITE(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++'
        END IF
        IF( PRESENT(INFO) ) THEN
          INFO = LINFO
        END IF
      END SUBROUTINE ERINFO

    
end module std_lapack
