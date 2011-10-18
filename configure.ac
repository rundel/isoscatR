AC_INIT([Rscat], 0.1.0)
AC_CONFIG_AUX_DIR([tools])

USE_MAGMA=yes

#################### GCC Compiler ####################

AC_PROG_CC(gcc)

#################### CUDA Toolkit ####################

AC_ARG_WITH(
	[cuda-home],
	[AC_HELP_STRING([--with-cuda-home=DIR],[full path to CUDA installation @<:@/usr/local/cuda@:>@])],
	[CUDA_HOME="$withval"],
	[CUDA_HOME=/usr/local/cuda]
)
AC_MSG_RESULT([setting CUDA home directory to ${CUDA_HOME}])

CUDA_INC=${CUDA_HOME}/include
AC_MSG_CHECKING([for CUDA headers])
AC_CHECK_FILE(
	[${CUDA_INC}/cublas.h],
	[CUDA_CPPFLAGS="-I${CUDA_INC}"],
	[USE_MAGMA=no]
)

AC_MSG_CHECKING([whether this is a 64-bit version of CUDA])
AC_CHECK_FILE(
	[${CUDA_HOME}/lib64/libcublas.so], 
	[CUDA_LDFLAGS="-L${CUDA_HOME}/lib64"],
	[CUDA_LDFLAGS="-L${CUDA_HOME}/lib"]
)

CUDA_LIBS="-lcuda -lcublas"

#################### MAGMA Libraries ####################

AC_ARG_WITH([magma-home],
	[AC_HELP_STRING([--with-magma-home=DIR],[full path to MAGMA installion directory @<:@/usr/local/magma@:>@])],
	[MAGMA_HOME="$withval"],
	[MAGMA_HOME=/usr/local/magma]
)

MAGMA_LIB=${MAGMA_HOME}/lib
MAGMA_INC=${MAGMA_HOME}/include

AC_MSG_CHECKING([for MAGMA headers])
AC_CHECK_FILE(
	[${MAGMA_INC}/magma.h],
	[],
	[USE_MAGMA=no]
)

AC_MSG_CHECKING([for MAGMA shared libraries])
AC_CHECK_FILE(
	[${MAGMA_LIB}/libmagma.so],
	[],
	[USE_MAGMA=no]
)

if test "${USE_MAGMA}" == "yes"; then
	CPPFLAGS="${CPPFLAGS} ${CUDA_CPPFLAGS} -DUSEMAGMA -I${MAGMA_INC}"
	LDFLAGS="${LDFLAGS} ${CUDA_LDFLAGS} -L${MAGMA_LIB}"
	LIBS="${LIBS} ${CUDA_LIBS} -lmagmablas -lmagma"
else
	AC_MSG_RESULT([Unable to locate magma library / headers, compiling without gpu support])
fi

#################### Makevars ####################

AC_SUBST([CPPFLAGS])
AC_SUBST([LDFLAGS])
AC_SUBST([LIBS])

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
