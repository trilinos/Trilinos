dnl @synopsis TAC_CHECK_FILE( file, [action-if-found], [action-if-not-found])
dnl
dnl identical to AC_CHECK_FILE except when cross-compiling.
dnl When cross-compiling, AC_CHECK_FILE aborts, 
dnl When cross-compiling TAC_CHECK_FILE does nothing
dnl
dnl @author Ken Stanley <ksstanl@sandia.gov>
dnl
AC_DEFUN([TAC_CHECK_FILE],
[
  if test "X$cross_compiling" != "Xyes"; then
    { AC_CHECK_FILE([$1],[$2],[$3]) } 
  fi
])


