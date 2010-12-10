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
zoltan_id_type="uint"
AC_ARG_WITH(id-type,
AC_HELP_STRING([--with-id-type], [Zoltan global ID type: unsigned int (default), int, long, or llong]),
[
if test "X$withval" == "Xint" ; then
  AC_DEFINE([HAVE_ZOLTAN_ID_TYPE_INT], [1], [Define if "typedef int ZOLTAN_ID_TYPE"])
  zoltan_id_type="int"
else
  if test "X$withval" == "Xllong" ; then
    AC_DEFINE([HAVE_ZOLTAN_ID_TYPE_LONG_LONG], [1], [Define if "typedef long long ZOLTAN_ID_TYPE"])
    zoltan_id_type="long long"
  else
    if test "X$withval" == "Xlong" ; then
      AC_DEFINE([HAVE_ZOLTAN_ID_TYPE_LONG], [1], [Define if "typedef long ZOLTAN_ID_TYPE"])
      zoltan_id_type="long"
    else
      if test "X$withval" == "Xuint" ; then
        AC_DEFINE([HAVE_ZOLTAN_ID_TYPE_UINT], [1], [Define if "typedef unsigned int ZOLTAN_ID_TYPE"])
      else
        AC_MSG_ERROR([Valid global ID types for Zoltan are uint, int, long, and llong])
      fi
    fi
  fi
fi
],
[
AC_DEFINE([HAVE_ZOLTAN_ID_TYPE_UINT], [1], [Define if "typedef unsigned int ZOLTAN_ID_TYPE"])
]
)
AC_MSG_RESULT([typedef $zoltan_id_type ZOLTAN_ID_TYPE])
]
)
