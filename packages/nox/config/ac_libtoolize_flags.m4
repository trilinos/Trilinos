dnl @synopsis AC_LIBTOOLIZE_COMPILER_FLAGS(COMPILER-FLAGS-VAR)
dnl
dnl Change the contents of variable COMPILER-FLAGS-VAR so that they are
dnl Libtool friendly, ie. prefix each of them with `-Xcompiler' so that
dnl Libtool doesn't remove them.
dnl
dnl @category Misc
dnl @author Ludovic Courtès <ludo@chbouib.org>
dnl @version 2004-09-07
dnl @license AllPermissive

AC_DEFUN([AC_LIBTOOLIZE_COMPILER_FLAGS],
  [ac_libtoolize_flags_temp=""
   for i in $$1
   do
     ac_libtoolize_flags_temp="$ac_libtoolize_flags_temp -Xcompiler $i"
   done
   $1="$ac_libtoolize_flags_temp"])dnl

dnl @synopsis AC_LIBTOOLIZE_LINKER_FLAGS(LINKER-FLAGS-VAR)
dnl
dnl Change the contents of variable LINKER-FLAGS-VAR so that they are
dnl Libtool friendly, ie. prefix each of them with `-Xlinker' so that
dnl Libtool doesn't remove them.
dnl
dnl @category Misc
dnl @author Ludovic Courtès <ludo@chbouib.org>
dnl @version 2004-09-07
dnl @license AllPermissive

AC_DEFUN([AC_LIBTOOLIZE_LINKER_FLAGS],
  [ac_libtoolize_flags_temp=""
   for i in $$1
   do
     ac_libtoolize_flags_temp="$ac_libtoolize_flags_temp -Xlinker $i"
   done
   $1="$ac_libtoolize_flags_temp"])dnl

dnl @synopsis AC_LIBTOOLIZE_CCLINKER_FLAGS(CCLINKER-FLAGS-VAR)
dnl
dnl Change the contents of variable CCLINKER-FLAGS-VAR so that they are
dnl Libtool friendly, ie. prefix each of them with `-XCClinker' so that
dnl Libtool doesn't remove them.
dnl
dnl @category Misc
dnl @author Ludovic Courtès <ludo@chbouib.org>
dnl @version 2004-09-07
dnl @license AllPermissive

AC_DEFUN([AC_LIBTOOLIZE_CCLINKER_FLAGS],
  [ac_libtoolize_flags_temp=""
   for i in $$1
   do
     ac_libtoolize_flags_temp="$ac_libtoolize_flags_temp -XCClinker $i"
   done
   $1="$ac_libtoolize_flags_temp"])dnl
