dnl @synopsis TAC_ARG_WITH_LIBS
dnl
dnl Test for --with-libs="name(s)".
dnl 
dnl Prepends the specified name(s) to the list of libraries to link 
dnl with.  
dnl
dnl Example use
dnl
dnl TAC_ARG_WITH_LIBS
dnl 
dnl tests for --with-libs and pre-pends to LIBS
dnl
dnl @author Jim Willenbring <jmwille@sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_LIBS],
[
AC_MSG_CHECKING([whether additional libraries are needed])
AC_ARG_WITH(libs,
AC_HELP_STRING([--with-libs], 
[List additional libraries here.  For example, --with-libs=-lsuperlu
or --with-libs=/path/libsuperlu.a]),
[
LIBS="${withval} ${LIBS}"
AC_MSG_RESULT([LIBS = ${LIBS}])
],
AC_MSG_RESULT(no)
)
]
)
