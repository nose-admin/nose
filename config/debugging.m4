AC_DEFUN([NOSE_SET_DEBUGGING],
	[

#
# look for debugging option
#
acx_debugging_ok=disable

AC_ARG_ENABLE(debugging,
        [AC_HELP_STRING([--enable-debugging=<lib>], [debugging <lib>])])
case $enable_debugging in
        yes ) acx_debugging_ok=yes ;;
        no | "") acx_debugging_ok=disable ;;
	*)  ;;
esac

if test $acx_debugging_ok = "yes"; then
	DEBUGGING=yes
	AC_SUBST(DEBUGGING)
else
	DEBUGGING=no
	AC_SUBST(DEBUGGING)
fi



#AC_DEFUN([NOSE_WITH_DEBUGGING],
#	[AC_ARG_WITH(debugging,
#	[ --with-debugging        compiles with debugging options],
#	[
#		echo "compilation mode set to debugging"	
#	])
#	AM_CONDITIONAL(WITH_DEBUGGING, test x"${with_debugging-no}" != xno)
#	]
#)

]
)