dnl @synopsis ZAC_ARG_WITH_ID
dnl
dnl Test for "--with-id-type="
dnl Default is "unsigned int".  Can also be "long", "llong" or "int".
dnl 
dnl Generates config.h macro.
dnl
AC_DEFUN([ZAC_ARG_WITH_ID],
[
AC_MSG_CHECKING([data type for ZOLTAN_ID_TYPE])
zoltan_id_type="unset"
AC_ARG_WITH(id-type,
AC_HELP_STRING([--with-id-type], [Zoltan global ID type: uint (default), ulong, or ullong]),
[
if test "X$withval" == "Xuint" ; then
  AC_DEFINE([UNSIGNED_INT_GLOBAL_IDS],[1],[define if ZOLTAN_ID_TYPE is unsigned int])
  zoltan_id_type="unsigned int"
else
  if test "X$withval" == "Xulong" ; then
    AC_DEFINE([UNSIGNED_LONG_GLOBAL_IDS],[1],[define if ZOLTAN_ID_TYPE is unsigned long])
    zoltan_id_type="unsigned long"
  else
    if test "X$withval" == "Xullong" ; then
      AC_DEFINE([UNSIGNED_LONG_LONG_GLOBAL_IDS],[1],[define if ZOLTAN_ID_TYPE is unsigned long long])
      zoltan_id_type="unsigned long long"
    else
       AC_MSG_ERROR([Valid global ID types for Zoltan are uint, ulong, and ullong])
    fi
  fi
fi
],
[
AC_DEFINE([UNSIGNED_INT_GLOBAL_IDS],[1],[define if ZOLTAN_ID_TYPE is unsigned int])
zoltan_id_type="unsigned int"
]
)
AC_MSG_RESULT([typedef $zoltan_id_type ZOLTAN_ID_TYPE])
]
)
