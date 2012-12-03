#
# Fortran flags setting
#
#
#
AC_DEFUN([NOSE_SET_FORTRAN_FLAGS],
	[

AC_REQUIRE([ACX_MPI])
AC_REQUIRE([NOSE_SET_DEBUGGING])
AC_REQUIRE([NOSE_SET_CUDA_FLAGS])

#FC=gfortran
#FC=ifort

#
# Set the flags according to the platform and compiler
#
case $host in

	x86*)
		case ${FC} in
			f95|g95)
				F77=${FC}
				FCFLAGS=" -ffloat-store -msse2 -O3"
				case ${DEBUGGING} in
					yes)
						FCFLAGS="-g"
						;;
					no)
						FCFLAGS=${FCFLAGS}
						;;
					*);;
				esac
				;;
			gfortran)
				F77="gfortran"	
				FC=$F77	
				FCFLAGS=" -ffloat-store  -ffree-line-length-none -msse2 -O3"
				case ${DEBUGGING} in
					yes)
						FCFLAGS="-g -fbounds-check -ffree-line-length-none"
						;;
					no)
						FCFLAGS=${FCFLAGS}
						;;
					*);;
				esac
				;;
			pgf77)
				F77="pgf95"
				FC=$F77	
				FCFLAGS="-O3  "
				case ${DEBUGGING} in
					yes)
						FCFLAGS="-g -ffree-line-length-none"
						;;
					no)
						FCFLAGS=${FCFLAGS}
						;;
					*);;
				esac
				;;
			ifort)
				F77="ifort"	
				FC=$F77	
				FCFLAGS=" -O3"
				case ${DEBUGGING} in
					yes)
						FCFLAGS="-g "
						;;
					no)
						FCFLAGS=${FCFLAGS}
						;;
					*);;
				esac
				;;
			*);;
		esac;;
		
    i?86*)    
 		case ${FC} in
			f95|g95)
				F77=${FC}
				FCFLAGS="-ffloat-store -ffree-line-length-none -msse2 -O3 -cpp -fsloppy-char"
				case ${DEBUGGING} in
					yes)
						FCFLAGS="-g"
						;;
					no)
						FCFLAGS=${FCFLAGS}
						;;
					*);;
				esac
				;;
			gfortran)
				F77="gfortran"
				FC=$F77	
				FCFLAGS="-ffloat-store -O3 -ffree-line-length-none "
				case ${DEBUGGING} in
					yes)
						FCFLAGS="-g -ffree-line-length-none"
						;;
					no)
						FCFLAGS=${FCFLAGS}
						;;
					*);;
				esac
				;;
			*);;
			pgf77)
				F77="pgf95"
				FC=$F77	
				FCFLAGS="-O3  "
				case ${DEBUGGING} in
					yes)
						FCFLAGS="-g -ffree-line-length-none"
						;;
					no)
						FCFLAGS=${FCFLAGS}
						;;
					*);;
				esac
				;;
			*);;
		esac;;
    	    	
esac

#
# if MPI not present, set serial compilation
#

if test x = x${MPIFC}; then
	MPIFC=${F77}
fi

if test ${USEMPI} = yes; then
	FCFLAGS="${FCFLAGS} -DHAVE_MPI"
fi

AC_SUBST([FCFLAGS])
AC_SUBST([FCFLAGS_F])
AC_SUBST([FCFLAGS_F90])
AC_SUBST([F77])
AC_SUBST([MPIFC])
			
	]
)
