dnl @synopsis TAC_ARG_CONFIG_EXPAT
dnl
dnl Support for the expat XML parser, optionally used by TSF.
dnl
dnl Test a variety of expat options:
dnl --enable-expat       - Turns expat XML parsing
dnl --with-expat         - specify root directory of expat
dnl --with-expat-include - specify include directory for expat 
dnl --with-expat-libs    - specify expat libraries
dnl --with-expat-libdir  - specify location of expat libraries
dnl
dnl If any of these options are set, HAVE_EXPAT will be defined for both
dnl Autoconf and Automake, and HAVE_EXPAT will be defined in the
dnl generated config.h file
dnl
dnl Shamelessly copied from Mike Heroux's macros for configuring MPI
dnl
dnl @author Kevin Long (krlong@sandia.gov)
dnl
AC_DEFUN([TAC_ARG_CONFIG_EXPAT],
[
AC_ARG_ENABLE(expat,
AC_HELP_STRING([--enable-expat],[Enable expat XML parser]),
[
HAVE_PKG_EXPAT=yes
EXPAT_LIBS=-lexpat
],
[HAVE_PKG_EXPAT=no]
)

AC_ARG_WITH(expat,
AC_HELP_STRING([--with-expat],[specify root directory of EXPAT installation]),
[
HAVE_PKG_EXPAT=yes
EXPAT_DIR=${withval}
AC_MSG_CHECKING(EXPAT directory)
AC_MSG_RESULT([${EXPAT_DIR}])
]
)

AC_ARG_WITH(expat-include,
AC_HELP_STRING([--with-expat-include],[specify include directory for EXPAT]),
[
HAVE_PKG_EXPAT=yes
EXPAT_INC=${withval}
AC_MSG_CHECKING(user-defined EXPAT includes)
AC_MSG_RESULT([${EXPAT_INC}])
]
)

AC_ARG_WITH(expat-libs,
AC_HELP_STRING([--with-expat-libs],[specify EXPAT libraries]),
[
HAVE_PKG_EXPAT=yes
EXPAT_LIBS=${withval}
AC_MSG_CHECKING(user-defined EXPAT libraries)
AC_MSG_RESULT([${EXPAT_LIBS}])
]
)

AC_ARG_WITH(expat-libdir,
AC_HELP_STRING([--with-expat-libdir],[specify location of EXPAT libraries]),
[
HAVE_PKG_EXPAT=yes
EXPAT_LIBDIR=${withval}
AC_MSG_CHECKING(user-defined EXPAT libraries)
AC_MSG_RESULT([${EXPAT_LIBS}])
]
)

AC_MSG_CHECKING(whether we are using EXPAT)
AC_MSG_RESULT([${HAVE_PKG_EXPAT}])

if test "X${HAVE_PKG_EXPAT}" = "Xyes"; then
   AC_DEFINE(HAVE_EXPAT,,[define if we want to use EXPAT])
fi

AC_SUBST(EXPAT_LIBS)
AC_SUBST(EXPAT_INC)


dnl Define Automake version of HAVE_EXPAT if appropriate

AM_CONDITIONAL(HAVE_EXPAT, [test "X${HAVE_PKG_EXPAT}" = "Xyes"])
])
