dnl @synopsis TAC_ARG_WITH_LAPACKLIB
dnl
dnl Test for --with-lapacklib="name".
dnl 
dnl Prepends the specified name to the list of files to check for LAPACK
dnl routines.  
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_LAPACKLIB],
[
AC_ARG_WITH(lapacklib,
AC_HELP_STRING([--with-lapacklib], 
[name of library containing LAPACK: will search lib directories for -lname]),
[
USE_LAPACKLIB=yes
NEWLAPACKLIB=${withval}
]
)

   LAPACKLIBS="cxml lapack complib.sgimath"

if test "X${USE_LAPACKLIB}" = "Xyes"; then

   LAPACKLIBS="${NEWLAPACKLIB} ${LAPACKLIBS}"
fi
])

