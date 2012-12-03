#
# Macros to support CUDA technology
#
#
# Last modified 2009/02/26
#
# Tomas Mancal,  mancal@karlov.mff.cuni.cz
#
#
AC_DEFUN([NOSE_SET_CUDA_FLAGS],
	[

AC_REQUIRE([NOSE_SET_DEBUGGING])
#
# This piece of code (starting at YvainCode -> and finishing at <- YvainCode)
# has be havily modified, but is based on a NVIDIA forum post by user Yvain
#
# YvainCode ->
# ------------------------------------------------------------------------------
# Setup CUDA paths
# ------------------------------------------------------------------------------
AC_ARG_WITH([cuda],
   [  --with-cuda=PATH    prefix where cuda is installed [default=auto]])
if test -n "$with_cuda"
then
   CUDA_CFLAGS="-I$with_cuda/include"
   CUDA_LIBS="-L$with_cuda/lib -lcudart"
   NVCC="$with_cuda/bin/nvcc"
   USE_CUDA=true
   AC_SUBST(CUDA_CFLAGS)
   AC_SUBST(CUDA_LIBS)
   AC_SUBST(NVCC)
   AC_DEFINE(USE_CUDA,1,[Define if we use CUDA.])   
else
#   CUDA_CFLAGS="-I/usr/local/cuda/include"
#   CUDA_LIBS="-L/usr/local/cuda/lib -lcuda -lcudart"
#   NVCC="nvcc"
   USE_CUDA=false
fi

AM_CONDITIONAL(CUDA, test x$USE_CUDA = xtrue )
AC_SUBST(USE_CUDA)

AC_ARG_ENABLE([emu],
   [  --enable-emu    Turn on device emulation for CUDA],
   [case "${enableval}" in
	yes) EMULATION=true;;
	no)  EMULATION=false;;
	*) AC_MSG_ERROR([bad value ${enableval} for --enable-emu]);;
   esac],
   [EMULATION=false]
)

# ------------------------------------------------------------------------------
# Setup nvcc flags
# ------------------------------------------------------------------------------
if test x$USE_CUDA = xtrue
then

	if test x$DEBUGGING = xyes
	then
   		NVCCFLAGS="-g"
	else
   		NVCCFLAGS="-O3 -use_fast_math"
	fi
	if test x$EMULATION = xtrue
	then
   		NVCCFLAGS+=" -deviceemu"
	fi

	AC_SUBST(NVCCFLAGS)

fi
# <-YvainCode

]

)