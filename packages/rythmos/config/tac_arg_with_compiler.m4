dnl @synopsis TAC_ARG_WITH_COMPILER(lcase_name, UCASE_NAME)
dnl
dnl Test for --with-lcase_name="compiler".  if defined, replace 
dnl the current UCASE_NAME definition.
dnl
dnl Use this macro to specify a compiler rather than letting 
dnl a built-in macro choose.
dnl
dnl NOTE: Only use this option for the serial case - if compiling
dnl with '--enable-mpi', use '--with-mpi-lcase_name="compiler"
dnl instead.
dnl
dnl Example use
dnl 
dnl TAC_ARG_WITH_COMPILER(cxx, CXX)
dnl 
dnl tests for --with-cxx and, if found, relpaces the existing
dnl value of CXX
dnl
dnl @author Jim Willenbring <jmwille@sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_COMPILER],
[
AC_MSG_CHECKING([whether a user defined [$2] compiler should be used])
AC_ARG_WITH($1,
AC_HELP_STRING([--with-$1], 
[user defined [$2] compiler found: will replace existing value of [$2]]),
[
$2="${withval}"
AC_MSG_RESULT([$2 = ${$2}])
],
AC_MSG_RESULT(no)
)
])

