dnl @synopsis TAC_ARG_CONFIG_MPI
dnl
dnl Test a variety of MPI options:
dnl --enable-mpi       - Turns MPI compiling mode on
dnl --with-mpi         - specify root directory of MPI
dnl --with-mpi-cxx     - specify MPI C++ compiler
dnl --with-mpi-cc      - specify MPI C compiler
dnl --with-mpi-f77     - specify MPI Fortran 77 compiler
dnl --with-mpi-include - specify include directory for MPI 
dnl --with-mpi-libs    - specify MPI libraries
dnl --with-mpi-libdir  - specify location of MPI libraries
dnl
dnl If any of these options are set, HAVE_MPI will be defined for both
dnl Autoconf and Automake, and HAVE_MPI will be defined in the
dnl generated config.h file
dnl
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_CONFIG_MPI],
[
AC_ARG_ENABLE(mpi,
AC_HELP_STRING([--enable-mpi],[enable MPI using mpiCC C++ compiler]),
[
HAVE_PKG_MPI=yes
MPI_CXX=mpiCC
MPI_CC=mpicc
MPI_F77=mpif77
],
[HAVE_PKG_MPI=no]
)

AC_ARG_WITH(mpi,
AC_HELP_STRING([--with-mpi],[specify root directory of MPI installation]),
[
HAVE_PKG_MPI=yes
MPI_DIR=${withval}
AC_MSG_CHECKING(MPI directory)
AC_MSG_RESULT([${MPI_DIR}])
]
)

AC_ARG_WITH(mpi-cxx,
AC_HELP_STRING([--with-mpi-cxx],[specify MPI C++ compiler]),
[
HAVE_PKG_MPI=yes 
MPI_CXX=${withval}
AC_MSG_CHECKING(user-defined MPI C++ compiler)
AC_MSG_RESULT([${MPI_CXX}])
]
)

AC_ARG_WITH(mpi-cc,
AC_HELP_STRING([--with-mpi-cc],[specify MPI C compiler]),
[
HAVE_PKG_MPI=yes 
MPI_CC=${withval}
AC_MSG_CHECKING(user-defined MPI C compiler)
AC_MSG_RESULT([${MPI_CC}])
]
)

AC_ARG_WITH(mpi-f77,
AC_HELP_STRING([--with-mpi-f77],[specify MPI Fortran 77 compiler]),
[
HAVE_PKG_MPI=yes 
MPI_F77=${withval}
AC_MSG_CHECKING(user-defined MPI Fortran 77 compiler)
AC_MSG_RESULT([${MPI_F77}])
]
)

AC_ARG_WITH(mpi-include,
AC_HELP_STRING([--with-mpi-include],[specify include directory for MPI]),
[
HAVE_PKG_MPI=yes
MPI_INC=${withval}
AC_MSG_CHECKING(user-defined MPI includes)
AC_MSG_RESULT([${MPI_INC}])
]
)

AC_ARG_WITH(mpi-libs,
AC_HELP_STRING([--with-mpi-libs],[specify MPI libraries]),
[
HAVE_PKG_MPI=yes
MPI_LIBS=${withval}
AC_MSG_CHECKING(user-defined MPI libraries)
AC_MSG_RESULT([${MPI_LIBS}])
]
)

AC_ARG_WITH(mpi-libdir,
AC_HELP_STRING([--with-mpi-libdir],[specify location of MPI libraries]),
[
HAVE_PKG_MPI=yes
MPI_LIBDIR=${withval}
AC_MSG_CHECKING(user-defined MPI libraries)
AC_MSG_RESULT([${MPI_LIBS}])
]
)

AC_MSG_CHECKING(whether we are using MPI)
AC_MSG_RESULT([${HAVE_PKG_MPI}])

if test "X${HAVE_PKG_MPI}" = "Xyes"; then
   AC_DEFINE(HAVE_MPI,,[define if we want to use MPI])
fi

dnl Define Automake version of HAVE_MPI if appropriate

AM_CONDITIONAL(HAVE_MPI, [test "X${HAVE_PKG_MPI}" = "Xyes"])


dnl
dnl --------------------------------------------------------------------
dnl Check for MPI compilers (must be done *before* AC_PROG_CXX,
dnl AC_PROG_CC and AC_PROG_F77)
dnl 
dnl --------------------------------------------------------------------

if test -n "${MPI_CXX}"; then
  AC_CHECK_PROG(MPI_CXX_EXISTS, ${MPI_CXX}, yes, no)
  if test "X${MPI_CXX_EXISTS}" = "Xyes"; then
    CXX=${MPI_CXX}
  else
    AC_MSG_ERROR([MPI C++ compiler (${MPI_CXX}) not found.])
  fi
fi

if test -n "${MPI_CC}"; then
  AC_CHECK_PROG(MPI_CC_EXISTS, ${MPI_CC}, yes, no)
  if test "X${MPI_CC_EXISTS}" = "Xyes"; then
    CC=${MPI_CC}
  else
    AC_MSG_ERROR([MPI C compiler (${MPI_CC}) not found.])
  fi
fi

if test -n "${MPI_F77}"; then
  AC_CHECK_PROG(MPI_F77_EXISTS, ${MPI_F77}, yes, no)
  if test "X${MPI_F77_EXISTS}" = "Xyes"; then
    F77=${MPI_F77}
  else
    AC_MSG_ERROR([MPI Fortran 77 compiler (${MPI_F77}) not found.])
  fi
fi
])
