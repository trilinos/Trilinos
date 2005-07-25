# TRILINOS_HEADER now takes three arguments:
# TRILINOS_HEADER(<library name>,<file name>,<subdirectory>)
# This is needed for multiple libraries being supplied by the
# same trilinos package (anasazi).
AC_DEFUN([TRILINOS_HEADER],
[
  save_CPPFLAGS=$CPPFLAGS
  done=no

  # First try the directories we already have
  CPPFLAGS="$save_CPPFLAGS $ABS_TRILINOS_CPPFLAGS"
  AC_PREPROC_IFELSE([AC_LANG_SOURCE(
  [[
#undef PACKAGE
#undef PACKAGE_NAME
#undef PACKAGE_BUGREPORT
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION
#include "$2"
  ]])],[done=yes],[])
  CPPFLAGS=$save_CPPFLAGS

  # Next try the Trilinos build dir
  if test $done = no; then
    CPPFLAGS="-I../$1/$3 $save_CPPFLAGS $ABS_TRILINOS_CPPFLAGS"
    AC_PREPROC_IFELSE([AC_LANG_SOURCE(
  [[
#undef PACKAGE
#undef PACKAGE_NAME
#undef PACKAGE_BUGREPORT
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION
#include "$2"
  ]])],[ done=yes
      REL_TRILINOS_CPPFLAGS="-I\$(top_builddir)/../$1/$3 ${REL_TRILINOS_CPPFLAGS}"
      ABS_TRILINOS_CPPFLAGS="-I../$1/$3 ${ABS_TRILINOS_CPPFLAGS}"
    ],[])
    CPPFLAGS=$save_CPPFLAGS
  fi

  # Next try the Trilinos source dir
  if test $done = no && test -n "$srcdir"; then
    CPPFLAGS="-I$srcdir/../$1/$3 $save_CPPFLAGS $ABS_TRILINOS_CPPFLAGS"
    AC_PREPROC_IFELSE([AC_LANG_SOURCE(
  [[
#undef PACKAGE
#undef PACKAGE_NAME
#undef PACKAGE_BUGREPORT
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION
#include "$2"
  ]])],[ done=yes
      REL_TRILINOS_CPPFLAGS="-I\$(top_srcdir)/../$1/$3 ${REL_TRILINOS_CPPFLAGS}"
      ABS_TRILINOS_CPPFLAGS="-I$srcdir/../$1/$3 ${ABS_TRILINOS_CPPFLAGS}"
    ],[])
    CPPFLAGS=$save_CPPFLAGS
  fi

  if test $done = no; then
    echo "------"
    echo "Cannot preprocess the include file $2, part of the $1 package."
    echo "Try --with-incdir=\"-I<dir1> -I<dir2>\""
    echo "------"
    AC_MSG_ERROR([Cannot find $1 include directory containing $2])
  fi
])
