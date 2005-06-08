dnl @synopsis TAC_ARG_WITH_AR
dnl
dnl Test for --with-ar="ar_program ar_flags".
dnl Default is "ar cru"
dnl 
dnl Generates an Automake conditional USE_ALTERNATE_AR that can be tested.  
dnl Generates the user-specified archiver command in @ALTERNATE_AR@.
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_AR],
[
AC_ARG_WITH(ar,
AC_HELP_STRING([--with-ar], [override archiver command (default is "ar cru")]),
[
AC_MSG_CHECKING(user-defined archiver)
AC_MSG_RESULT([${withval}])
USE_ALTERNATE_AR=yes
ALTERNATE_AR="${withval}"
]
)

if test -n "${SPECIAL_AR}" && test "X${USE_ALTERNATE_AR}" != "Xyes";
then
  USE_ALTERNATE_AR=yes
  ALTERNATE_AR="${SPECIAL_AR}"
fi

# CHANGES for LIBTOOL 
# Note:  We now overwrite the AR and AR_FLAGS variables with the alternate
# archiver.  This is necessary for libtool to respect the alternate archiver.
# Note that the AR_FLAGS is specifically set to a space instead of being
# empty because libtool.m4 automatically sets AR_FLAGS=cru if AR_FLAGS has
# zero length.  We don't want this when using an alternate archiver.
#
# It may no longer be necessary to explicitly specify ALTERNATE_AR in each
# Makefile.am with this change.
AC_MSG_CHECKING(for special archiver command)
if test "X${USE_ALTERNATE_AR}" = "Xyes"; then
   AC_MSG_RESULT([${ALTERNATE_AR}])
   AM_CONDITIONAL(USE_ALTERNATE_AR, true)
   AR="${ALTERNATE_AR}"
   AR_FLAGS=" "
else
   AC_MSG_RESULT([none])
   AM_CONDITIONAL(USE_ALTERNATE_AR, false)
   AC_CHECK_TOOL(AR, ar, false)
   if test "X${AR}" = "Xfalse" ; then
      AC_MSG_ERROR([Default archiver ar not found.  Specify an archiver with --with-ar configure flag.])
   fi
fi
AC_SUBST(ALTERNATE_AR)
])

