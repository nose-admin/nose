#
#
#
AC_INIT(nose,0.5.10.920,tomas.mancal@mff.cuni.cz)

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR([src])
AC_LANG(Fortran)

AC_CANONICAL_HOST

#ACX_F77_CMAIN_FFLAGS
AC_PROG_CC
AC_PROG_RANLIB
AC_PROG_LEX
AC_PROG_YACC
AX_PROG_TCL
AC_LANG(C)
ACX_BLAS
ACX_LAPACK
ACX_MPI
AC_LANG(Fortran)
ACX_MPI


NOSE_SET_DEBUGGING

NOSE_SET_CUDA_FLAGS

NOSE_SET_FORTRAN_FLAGS


AC_CONFIG_FILES([Makefile lib/Makefile lib/cuda/Makefile src/Makefile src/scripts/Makefile \
  src/pert/Makefile src/utils/Makefile src/pert/modules/Makefile \
  src/pert/modules/tdpt3/Makefile src/pert/modules/qme/Makefile src/pert/modules/test/Makefile \
  src/pert/modules/montecarlo/Makefile \
  src/pert/modules/molc/Makefile \
  testsuite/Makefile testsuite/test01/Makefile testsuite/test02/Makefile testsuite/test03/Makefile \
  testsuite/test04/Makefile testsuite/test05/Makefile \
  doc/Makefile ])

AM_INIT_AUTOMAKE
AC_OUTPUT


