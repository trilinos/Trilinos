dnl @synopsis ACX_GGEV
dnl
dnl This macro looks for the LAPACK generalized eigenvalue routine
dnl _GGEV (see http://www.netlib.org/lapack/).  There is a deprecated
dnl name for the same function, _GEGV, that is used instead on some
dnl architectures (namely DECs).  This macro attempts to determine which
dnl version is available and defines the following variables:
dnl HAVE_LAPACK_GGEV (true or false) and HAVE_LAPACK_GEGV (true or false).  
dnl This macro requires the LAPACK_LIBS, BLAS_LIBS, and FLIBS variables 
dnl set by the TAC_LAPACK macro.
dnl
dnl NOTE:  This is a modification of Eric Phipps' script for NOX.
dnl
dnl @author Heidi Thornquist <etphipp@sandia.gov>

AC_DEFUN([ACX_GGEV], [
acx_ggev_ok=no
acx_gegv_ok=no

save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

# This test can be performed for one data type (double precision) and generalized.
# First, try dggev
AC_F77_FUNC(dggev)
AC_MSG_CHECKING([for $dggev in $LAPACK_LIBS])
AC_TRY_LINK_FUNC($dggev, [acx_ggev_ok=yes])
AC_MSG_RESULT($acx_ggev_ok)
if test $acx_ggev_ok = yes; then
	AC_DEFINE(HAVE_LAPACK_GGEV,1,[Define if you have LAPACK _GGEV.])
else
	# Next try dgegv
	AC_F77_FUNC(dgegv)
	AC_MSG_CHECKING([for $dgegv in $LAPACK_LIBS])
	AC_TRY_LINK_FUNC($dgegv, [acx_gegv_ok=yes])
	AC_MSG_RESULT($acx_gegv_ok)
	if test $acx_gegv_ok = yes; then
		AC_DEFINE(HAVE_LAPACK_GEGV,1,[Define if you have LAPACK _GEGV.])	
	fi
fi
	
LIBS="$save_LIBS"

])dnl ACX_GGEV
