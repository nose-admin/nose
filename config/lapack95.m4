#-------------------------------------------------------------------------------------------------------
#
# Tries to guess location of lapack95 according to the platform
#
AC_DEFUN([NOSE_LAPACK95_GUESS_DIR],
	[	
	 	case $host in

		powerpc-ibm-*)
			LAPACK95_DIR="-L${PWD}/lib/ -L${PWD}/lib/lapack95"
			LAPACK95_INC=-I${PWD}/lib/lapack95/lapack95_modules/
			LAPACK95_LDADD="-llapack95"
			LAPACK95_LINK="make.inc_power2"
			;;

		x86*)
			LAPACK95_DIR="-L/usr/lib64/ -L${PWD}/lib/lapack95"
			LAPACK95_INC="-I${PWD}/lib/lapack95/lapack95_modules/"
			case ${FC} in
				f95|g95)
					LAPACK95_LDADD="-llapack95 -llapack_g95 -lblas -lg2c"
					LAPACK95_LINK="make.inc_linux_g95"
					;;
				gfortran)
					LAPACK95_LDADD="-llapack95"
					LAPACK95_LINK="make.inc_linux_gfortran"
					;;
				*);;
			esac;;
		
    	i?86*)    
    		LAPACK95_DIR="-L/usr/lib"
    		LAPACK95_INC="-I/usr/lib/lapack95_modules/"
 			case ${FC} in
				f95|g95)
					LAPACK95_LDADD="-llapack95 -llapack_g95 -lblas -lg2c"
					LAPACK95_LINK="make.inc_linux_g95"
					;;
				gfortran)
					LAPACK95_LDADD="-llapack95"
					LAPACK95_LINK="make.inc_linux_gfortran"
					;;
				*);;
			esac;;
    	
    	alpha*)
    		LAPACK95_DIR="-L${PWD}/lib/lapack95"
    		LAPACK95_INC="-I${PWD}/lib/lapack95/lapack95_modules/"
    		LAPACK95_LDADD="-llapack95 -ldxml"
			LAPACK95_LINK="make.inc_alpha"
    		;;
    	
	 esac	
])
	 
#-----------------------------------------------------------------------------------
#  Sets the lapack95 directory to the one provided by the package
#  and schedules its compilation
#
AC_DEFUN([NOSE_LAPACK95_SET_PACKAGE_DIR],
	[	
	 	case $host in

		powerpc-ibm-*)
			LAPACK95_DIR=-L${PWD}/lib/
			LAPACK95_INC=-I${PWD}/lib/lapack95_modules/
			LAPACK95_LDADD="-llapack95 -lessl"
			LAPACK95_LINK="make.inc_power2"
			;;

		x86*)
			LAPACK95_DIR="-L/usr/lib64/ -L${PWD}/lib/lapack95"
			LAPACK95_INC="-I${PWD}/lib/lapack95/lapack95_modules/"
			case ${FC} in
				f95|g95)
					LAPACK95_LDADD="-llapack95 -llapack_g95 -lblas -lg2c"
					LAPACK95_LINK="make.inc_linux_g95"
					;;
				gfortran)
					LAPACK95_LDADD="-llapack95"
					LAPACK95_LINK="make.inc_linux_gfortran"
					;;
				*);;
			esac;;
		
    	i?86*)    
    		LAPACK95_DIR="-L/usr/lib"
    		LAPACK95_INC="-I/usr/lib/lapack95_modules/"
 			case ${FC} in
				f95|g95)
					LAPACK95_LDADD="-llapack95 -llapack_g95 -lblas -lg2c"
					LAPACK95_LINK="make.inc_linux_g95"
					;;
				gfortran)
					LAPACK95_LDADD="-llapack95"
					LAPACK95_LINK="make.inc_linux_gfortran"
					;;
				*);;
			esac;;
    	
    	alpha*)
    		LAPACK95_DIR="-L${PWD}/lib/lapack95"
    		LAPACK95_INC="-I${PWD}/lib/lapack95/lapack95_modules/"
    		LAPACK95_LDADD="-llapack95 -ldxml"
			LAPACK95_LINK="make.inc_alpha"
    		;;
    	
	 esac	
])
	 

#---------------------------------------------------------------------------------------
#
# Checks if the lapack95 links with the current settings
#
AC_DEFUN([NOSE_CHECK_LINK_LAPACK95],
	[AC_REQUIRE([ACX_LAPACK])
	 AC_LINK_IFELSE(
	 	[AC_LANG_PROGRAM([],[[
	 		external func
	 		call func()
	 	]])],
	 	[
	 		AC_MSG_RESULT([yes])
	 	],
	 	[
	 		AC_MSG_RESULT([no])
	 		AC_MSG_RESULT([using lapack95 from this package])
	 	]
	 )]	 
)

#---------------------------------------------------------------------------------------
#
# Finds and sets lapack95
#
AC_DEFUN([NOSE_SUBST_LAPACK95],
	[	
	 	AC_LANG(Fortran)
	 	AC_MSG_CHECKING([for lapack95])
	 	NOSE_LAPACK95_GUESS_DIR
	 	NOSE_CHECK_LINK_LAPACK95
		AC_SUBST([LAPACK95_DIR])
		AC_SUBST([LAPACK95_INC])
		AC_SUBST([LAPACK95_LDADD])
		AC_SUBST([LAPACK95_LINK])
		AC_MSG_RESULT([default setting... LAPACK95_INC=]$LAPACK95_INC)
		AC_MSG_RESULT([                   LAPACK95_LINK=]$LAPACK95_LINK)
	]
)
