dnl @synopsis TAC_ARG_WITH_LIBDIRS
dnl
dnl Test for --with-libdirs="-Llibdir1 -Llibdir2".  if defined, 
dnl prepend "-Llibdir1 -Llibdir2" to LDFLAGS
dnl
dnl Use this macro to facilitate addition of directories to library search path.
dnl 
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_LIBDIRS],
[
AC_MSG_CHECKING([whether additional library search paths defined])
AC_ARG_WITH(libdirs,
AC_HELP_STRING([--with-libdirs], 
[OBSOLETE use --with-ldflags instead.  (ex. --with-ldflags="-L<DIR> -L<DIR2>")]),
[
LDFLAGS="${withval} ${LDFLAGS}"
AC_MSG_RESULT([${withval}])
],
AC_MSG_RESULT(no)
)
])

