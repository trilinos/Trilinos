dnl @synopsis AC_PROMPT_USER(VARIABLENAME,QUESTION,[DEFAULT])
dnl
dnl Asks a QUESTION and puts the results in VARIABLENAME with an optional
dnl DEFAULT value if the user merely hits return.  Also calls
dnl AC_DEFINE_UNQUOTED() on the VARIABLENAME for VARIABLENAMEs that should
dnl be entered into the config.h file as well.
dnl
dnl @version $Id$
dnl @author Wes Hardaker <wjhardaker@ucdavis.edu>
dnl
AC_DEFUN([AC_PROMPT_USER],
[
MSG_CHECK=`echo "$2" | tail -1`
AC_CACHE_CHECK($MSG_CHECK, ac_cv_user_prompt_$1,
[echo "" >&AC_FD_MSG
AC_PROMPT_USER_NO_DEFINE($1,[$2],$3)
eval ac_cv_user_prompt_$1=\$$1
echo $ac_n "setting $MSG_CHECK to...  $ac_c" >&AC_FD_MSG
])
if test "$ac_cv_user_prompt_$1" != "none"; then
  if test "$4" != ""; then
    AC_DEFINE_UNQUOTED($1,"$ac_cv_user_prompt_$1")
  else
    AC_DEFINE_UNQUOTED($1,$ac_cv_user_prompt_$1)
  fi
fi
]) dnl

