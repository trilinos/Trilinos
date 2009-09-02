dnl @synopsis ACX_DGGEV
dnl
dnl This macro looks for the LAPACK generalized eigenvalue routine
dnl DGGEV (see http://www.netlib.org/lapack/).  There is a deprecated
dnl name for the same function, DGEGV, that is used instead on some
dnl architectures (namely DECs).  This macro attempts to determine which
dnl version is available and defines the following variables:
dnl HAVE_LAPACK_GENEV (true or false), HAVE_LAPACK_DGGEV (true or false),
dnl and HAVE_LAPACK_DGEGV (true or false).  If neither function is found
dnl HAVE_LAPACK_GENEV is false.  This macro requires the LAPACK_LIBS, 
dnl BLAS_LIBS, and FLIBS variables set by the ACX_LAPACK macro.
dnl
dnl @version $Id$
dnl @author Eric Phipps <etphipp@sandia.gov>

AC_DEFUN([ACX_DGGEV], [
acx_dggev_ok=no
acx_dgegv_ok=no

save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

# First, try dggev
AC_F77_FUNC(dggev)
AC_MSG_CHECKING([for $dggev in $LAPACK_LIBS])
AC_TRY_LINK_FUNC($dggev, [acx_dggev_ok=yes])
AC_MSG_RESULT($acx_dggev_ok)
if test $acx_dggev_ok = yes; then
	AC_DEFINE(HAVE_LAPACK_GENEV,1,[Define if you have LAPACK DGGEV/DGEGV.])
	AC_DEFINE(HAVE_LAPACK_DGGEV,1,[Define if you have LAPACK DGGEV.])
else
	# Next try dgegv
	AC_F77_FUNC(dgegv)
	AC_MSG_CHECKING([for $dgegv in $LAPACK_LIBS])
	AC_TRY_LINK_FUNC($dgegv, [acx_dgegv_ok=yes])
	AC_MSG_RESULT($acx_dgegv_ok)
	if test $acx_dgegv_ok = yes; then
		AC_DEFINE(HAVE_LAPACK_GENEV,1,[Define if you have LAPACK DGGEV/DGEGV.])
		AC_DEFINE(HAVE_LAPACK_DGEGV,1,[Define if you have LAPACK DGEGV.])	
	fi
fi
	
LIBS="$save_LIBS"

])dnl ACX_DGGEV
