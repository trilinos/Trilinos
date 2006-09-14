dnl @synopsis AC_ARG_WITH_EXECINCLUDEDIR
dnl
dnl Test for --with-execincludedir="dir".  If defined, set 
dnl architecture-specific include directory to "dir", ${exec_prefix}/include
dnl otherwise.
dnl
dnl Use this macro to provide a path to where architecture-specific headers
dnl such as Sacado_config.h should go
dnl 
dnl
dnl @author Eric Phipps <etphipp@sandia.gov>
dnl
AC_DEFUN([AC_ARG_WITH_EXECINCLUDEDIR],
[
AC_MSG_CHECKING([for architecture-specific include install directory])
AC_ARG_WITH(execincludedir,
AC_HELP_STRING([--with-execincludedir], 
[path for installing architecture-specific header files]),
[
execincludedir=$withval
],
execincludedir='${exec_prefix}/include'
)
AC_MSG_RESULT([$execincludedir])
AC_SUBST([execincludedir])
])
