dnl @synopsis AC_PROMPT_USER_NO_DEFINE(VARIABLENAME,QUESTION,[DEFAULT])
dnl
dnl Asks a QUESTION and puts the results in VARIABLENAME with an optional
dnl DEFAULT value if the user merely hits return.
dnl
dnl @version $Id$
dnl @author Wes Hardaker <wjhardaker@ucdavis.edu>
dnl
AC_DEFUN([AC_PROMPT_USER_NO_DEFINE],
dnl changequote(<<, >>) dnl
dnl <<
[
if test "x$defaults" = "xno"; then
echo $ac_n "$2 ($3): $ac_c"
read tmpinput
if test "$tmpinput" = "" -a "$3" != ""; then
  tmpinput="$3"
fi
eval $1=\"$tmpinput\"
else
tmpinput="$3"
eval $1=\"$tmpinput\"
fi
]
dnl >>
dnl changequote([, ])
) dnl done AC_PROMPT_USER
