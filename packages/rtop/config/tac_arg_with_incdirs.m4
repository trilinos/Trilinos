dnl @synopsis TAC_ARG_WITH_INCDIRS
dnl
dnl Test for --with-incdirs="-Iincdir1 -Iincdir2".  if defined, prepend 
dnl "-Iincdir1 -Iincdir2" to CPPFLAGS
dnl
dnl Use this macro to facilitate addition of directories to include file search path.
dnl 
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_INCDIRS],
[
AC_MSG_CHECKING([whether additional include search paths defined])
AC_ARG_WITH(incdirs,
AC_HELP_STRING([--with-incdirs], 
[additional directories containing include files: will prepend to search here for includes, use -Idir format]),
[
CPPFLAGS="${withval} ${CPPFLAGS}"
AC_MSG_RESULT([${withval}])
],
AC_MSG_RESULT(no)
)
])

