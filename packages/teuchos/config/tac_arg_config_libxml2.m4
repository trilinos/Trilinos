dnl @synopsis TAC_ARG_CONFIG_LIBXML2
dnl
dnl Support for the libxml2 XML parser, optionally used by TSF.
dnl
dnl Test a variety of libxml2 options:
dnl --enable-libxml2       - Turns libxml2 XML parsing
dnl --with-libxml2         - specify root directory of libxml2
dnl --with-libxml2-include - specify include directory for libxml2 
dnl --with-libxml2-libdir  - specify location of libxml2 libraries
dnl
dnl If any of these options are set, HAVE_LIBXML2 will be defined for both
dnl Autoconf and Automake, and HAVE_LIBXML2 will be defined in the
dnl generated config.h file
dnl
dnl Shamelessly copied from Mike Heroux's macros for configuring MPI
dnl
dnl @author Kevin Long (krlong@sandia.gov)
dnl
AC_DEFUN([TAC_ARG_CONFIG_LIBXML2],
[
AC_ARG_ENABLE(libxml2,
AC_HELP_STRING([--enable-libxml2],[Enable libxml2 XML parser]),
[
HAVE_PKG_LIBXML2=yes
XML_LIBS=-lxml2
],
[HAVE_PKG_LIBXML2=no]
)

AC_ARG_WITH(libxml2-include,
AC_HELP_STRING([--with-libxml2-include],[specify include directory for LIBXML2]),
[
XML_INC=${withval}
AC_MSG_CHECKING(user-defined LIBXML2 includes)
AC_MSG_RESULT([${XML_INC}])
]
)

AC_ARG_WITH(libxml2-libdir,
AC_HELP_STRING([--with-libxml2-libdir],[specify location of LIBXML2 libraries]),
[
XML_LIBDIR=${withval}
AC_MSG_CHECKING(user-defined LIBXML2 libraries)
AC_MSG_RESULT([${XML_LIBDIRS}])
]
)

AC_MSG_CHECKING(whether we are using LIBXML2)
AC_MSG_RESULT([${HAVE_PKG_LIBXML2}])

if test "X${HAVE_PKG_LIBXML2}" = "Xyes"; then
   AC_DEFINE(HAVE_LIBXML2,,[define if we want to use LIBXML2])
fi

AC_SUBST(XML_LIBDIR)
AC_SUBST(XML_LIBS)
AC_SUBST(XML_INC)


dnl Define Automake version of HAVE_LIBXML2 if appropriate

AM_CONDITIONAL(HAVE_LIBXML2, [test "X${HAVE_PKG_LIBXML2}" = "Xyes"])
])
