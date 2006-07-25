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

AC_MSG_CHECKING(for special archiver command)
if test "X${USE_ALTERNATE_AR}" = "Xyes"; then
   AC_MSG_RESULT([${ALTERNATE_AR}])
   AM_CONDITIONAL(USE_ALTERNATE_AR, true)
else
   AC_MSG_RESULT([none])
   AM_CONDITIONAL(USE_ALTERNATE_AR, false)
fi
AC_SUBST(ALTERNATE_AR)
])

