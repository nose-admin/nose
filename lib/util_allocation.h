!
! Allocation counting macros
!

#define GNU_FORTRAN

#ifdef GNU_FORTRAN
#define A_SIZEOF(a) sizeof(a)
!!util_sizeof(kind(a))
#define B_SIZEOF(a) sizeof(a)
!!util_sizeof(kind(a),size(a))
#endif

!
! Allocation and deallocation of a scalar
!
#define ALLOCATE_S(a) call count_name("a"); allocate(a);call count_bytes(A_SIZEOF(a))  !;print *, get_mb_count()
#define DEALLOCATE_S(a) call uncount_name("a"); call uncount_bytes(A_SIZEOF(a)); deallocate(a)

!
! Allocation and deallocation of arrays
!
#define ALLOCATE(a,b) call count_name("a"); allocate(a b);call count_bytes(B_SIZEOF(a))
#define DEALLOCATE(a) call uncount_name("a"); call uncount_bytes(B_SIZEOF(a));deallocate(a)

!
! Uncounting of variables that were allocated inside a subroutine
! or function.
!
#define UNCOUNT(a) call uncount_bytes(B_SIZEOF(a)); call uncount_name("a")
