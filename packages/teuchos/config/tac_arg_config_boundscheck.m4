dnl @synopsis TAC_ARG_CONFIG_BOUNDSCHECK
dnl
dnl Support for boundschecking in Teuchos_Array
dnl
dnl --enable-boundscheck       - Turns on boundschecking
dnl
dnl @author Kevin Long (krlong@sandia.gov)
dnl
AC_DEFUN([TAC_ARG_CONFIG_BOUNDSCHECK],
[
AC_ARG_ENABLE(boundscheck,
AC_HELP_STRING([--enable-boundscheck],[enable Teuchos_Array boundschecking]),
[
BNDSCHK=yes
],
[BNDSCHK=no]
)

AC_MSG_CHECKING(whether we are using TSF array boundschecking)
AC_MSG_RESULT([${BNDSCHK}])

if test "X$BNDSCHK" = "Xyes"; then
  AC_DEFINE(TEUCHOS_HAVE_ARRAY_BOUNDSCHECK,,[Define if we are using array boundschecking])
fi

])

