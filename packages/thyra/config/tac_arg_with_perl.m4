dnl @synopsis TAC_ARG_WITH_PERL(DEFAULT_VAL)
dnl
dnl Test for --enable-gnumake and set to DEFAULT_VAL value if feature not specified.
dnl Calls AC_DEFINE to define HAVE_GNUMAKE if value is not equal to "no"
dnl Calls AM_CONDITIONAL to define USING_GNUMAKE to true/false.
dnl 
dnl This file was based on tac_arg_with_ar.m4 by Mike Heroux
dnl @author Roger Pawlowski <rppawlo@sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_PERL],
[

AC_ARG_WITH(perl,
AC_HELP_STRING([--with-perl], [supply a perl executable.  For example --with-perl=/usr/bin/perl.]),
[
AC_MSG_CHECKING(for user supplied perl executable)
AC_MSG_RESULT([${withval}])
USER_SPECIFIED_PERL=yes
PERL_EXE="${withval}"
],
[
USER_SPECIFIED_PERL=no
])

if test "X${USER_SPECIFIED_PERL}" = "Xyes"; then
  AC_CHECK_FILE(${PERL_EXE}, [HAVE_PERL=yes], [HAVE_PERL=no])
  AC_SUBST(PERL_EXE, ${PERL_EXE})
else
  AC_CHECK_PROG(HAVE_PERL, perl, yes, no)
  AC_SUBST(PERL_EXE, perl)
fi
AM_CONDITIONAL(USING_PERL, test X${HAVE_PERL} = Xyes)
])

