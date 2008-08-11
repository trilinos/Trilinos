dnl RYTHMOS Interface 
dnl
dnl @author Roger Pawlowski <rppawlo@sandia.gov>
dnl 09/13/06 tscoffe:  copied from nox
dnl
AC_DEFUN([RYTHMOS_INTERFACE],
[
  AC_ARG_ENABLE($1,
  [AC_HELP_STRING([--enable-$3-$1],[compile $3-$1 interface libraries])],
  [ADDON_$1=$enableval],
  [ADDON_$1=no]
  )

  AC_MSG_CHECKING(whether we should build the $3-$1 interface libraries)
  AC_MSG_RESULT(${ADDON_$1})
 
  AM_CONDITIONAL(BUILD_$2, test X${ADDON_$1} = Xyes)

])
