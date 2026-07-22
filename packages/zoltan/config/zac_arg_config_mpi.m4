dnl @synopsis ZAC_ARG_CONFIG_MPI
dnl
dnl Test a variety of MPI options:
dnl --enable-mpi       - Turns MPI compiling mode on
dnl --with-mpi         - specify root directory of MPI
dnl --with-mpi-compilers - Turns on MPI compiling mode and sets the MPI C++
dnl                       compiler C, and Fortran
dnl --with-mpi-incdir - specify include directory for MPI 
dnl --with-mpi-libs    - specify MPI libraries
dnl --with-mpi-libdir  - specify location of MPI libraries
dnl
dnl If any of these options are set, HAVE_MPI will be defined for both
dnl Autoconf and Automake, and HAVE_MPI will be defined in the
dnl generated config.h file
dnl
dnl if --disable-mpi, then Zoltan will build serial MPI.
dnl
dnl --enable-mpi and --with-mpi-compilers are the default.
dnl
dnl --without-mpi is actually a user error, but we'll interpret it as --disable-mpi
dnl
dnl If CC, CXX, F77 and/or F90/FTN/FC have been set by the user, and MPI compilers
dnl are desired, these will be assumed to be the MPI compilers.
dnl
dnl This was adapted from the Trilinos TAC_ARG_CONFIG_MPI.
dnl
AC_DEFUN([ZAC_ARG_CONFIG_MPI],
[

HAVE_PKG_MPI=unset
SEEK_MPI_COMPILERS=unset
MPI_COMPILER_PATH=unset

AC_ARG_ENABLE(mpi,
[AC_HELP_STRING([--enable-mpi],[enable MPI support])],
[
  if test X${enableval} = Xno; then
    HAVE_PKG_MPI=no
  else
    HAVE_PKG_MPI=yes
  fi
]
)

AC_ARG_WITH(mpi,
[AC_HELP_STRING([--with-mpi=MPIROOT],[the MPI root directory (above bin,lib,include), enables MPI])],
[
  if test X${withval} = Xno; then
    HAVE_PKG_MPI=no
  else
    HAVE_PKG_MPI=yes
    if test X${withval} != Xyes; then
      MPI_DIR=${withval}
    fi
  fi
]
)

AC_ARG_WITH(mpi-compilers,
[AC_HELP_STRING([--with-mpi-compilers={yes/no/path}],[Find MPI compilers/Don't use MPI compilers/Find MPI compilers in path])],
[
  HAVE_PKG_MPI=yes
  if test X${withval} = Xno; then
    SEEK_MPI_COMPILERS=no
  else
    SEEK_MPI_COMPILERS=yes
    if test X${withval} != Xyes; then
      MPI_COMPILER_PATH=${withval}
    fi
  fi
],
[
  if test X${HAVE_PKG_MPI} != Xno; then
    SEEK_MPI_COMPILERS=yes
    HAVE_PKG_MPI=yes
  fi
]
)

dnl Using MPI is the default

if test X${HAVE_PKG_MPI} = unset ; then
  HAVE_PKG_MPI=yes
fi

if test X${SEEK_MPI_COMPILERS} = Xyes; then

  if test X${MPI_COMPILER_PATH} != Xunset ; then
#    MPI_SEEK_PATH=$MPI_COMPILER_PATH$PATH_SEPARATOR$PATH
    MPI_SEEK_PATH=$MPI_COMPILER_PATH
  elif test -n "${MPI_DIR}" ; then
#    MPI_SEEK_PATH=$MPI_DIR/bin$PATH_SEPARATOR$PATH
    MPI_SEEK_PATH=$MPI_DIR/bin
  else
    MPI_SEEK_PATH=$PATH
  fi

  dnl Find C MPI compiler if MPI_CC is not already defined

  if test -z "${MPI_CC}"; then
    if test -f "${CC}"; then

      MPI_CC=${CC}

    else

      if test -n "${CC}" ; then
        MPI_CC_CANDIDATE=${CC}
      else
        MPI_CC_CANDIDATE=mpicc
      fi

      AC_PATH_PROG(MPI_CC, ${MPI_CC_CANDIDATE}, [notFound], [PATH = ${MPI_SEEK_PATH}])
  
      if test "${MPI_CC}" != "notFound" ; then
        CC=${MPI_CC}
      else
        echo "-----"
        echo "Cannot find MPI C compiler in " ${MPI_SEEK_PATH}
        echo "Specify a path to all mpi compilers with --with-mpi-compilers=PATH,"
        echo "or specify a path to top mpi directory (above bin) with --with-mpi=PATH,"
        echo "or specify a C compiler using CC=<compiler>"
        echo "or --disable-mpi"
        echo "-----"
        AC_MSG_ERROR([MPI C compiler not found.])
      fi
    fi
  fi

  if test "X$ac_cv_use_zoltan_cppdriver" = "Xyes"; then
    dnl Find C++ MPI compiler if MPI_CXX is not already defined
  
    if test -z "${MPI_CXX}"; then

      if test -f "${CXX}"; then
  
        MPI_CXX=${CXX}
  
      else
  
        if test -n "${CXX}" ; then
          MPI_CXX_CANDIDATES=${CXX}
        else
          MPI_CXX_CANDIDATES="[mpicxx mpic++ mpiCC]"
        fi
  
        AC_PATH_PROGS(MPI_CXX, ${MPI_CXX_CANDIDATES}, [notFound], [PATH = ${MPI_SEEK_PATH}])
    
        if test "${MPI_CXX}" != "notFound" ; then
          CXX=${MPI_CXX}
        else
          echo "-----"
          echo "Cannot find MPI C++ compiler in " ${MPI_SEEK_PATH}
          echo "Specify a path to all mpi compilers with --with-mpi-compilers=PATH,"
          echo "or specify a path to top mpi directory (above bin) with --with-mpi=PATH,"
          echo "or specify a C++ compiler using CXX=<compiler>"
          echo "or --disable-mpi"
          echo "-----"
          AC_MSG_ERROR([MPI C++ compiler not found.])
        fi
      fi
    fi
  fi

#  if test "X$ac_cv_use_fortran" = "Xyes"; then
#    dnl Find a Fortran 77 MPI compiler if MPI_F77 is not already defined
#  
#    if test -z "${MPI_F77}"; then
#      MPI_F77_CANDIDATE=mpif77
#      if test -n "${F77}"; then
#        MPI_F77_CANDIDATE=${F77}
#      fi
#  
#      AC_PATH_PROG(MPI_F77, ${MPI_F77_CANDIDATE}, [notFound], [PATH = ${MPI_SEEK_PATH}])
#  
#      if test "${MPI_F77}" != "notFound" ; then
#        F77=${MPI_F77}
#      else
#        echo "-----"
#        echo "Cannot find MPI Fortran 77 compiler."
#        echo "Specify a path to all mpi compilers with --with-mpi-compilers=PATH,"
#        echo "or specify a path to top mpi directory (above bin) with --with-mpi=PATH,"
#        echo "or specify a fortran 77 compiler using F77=<compiler>"
#        echo "-----"
#        AC_MSG_ERROR([MPI C compiler not found.])
#      fi
#    fi
#  fi

  if test "X$ac_cv_use_fortran90" = "Xyes"; then
    dnl Find a Fortran 90 MPI compiler if MPI_F90 is not already defined
  
    if test -z "${MPI_FC}"; then

      if test -f "${FC}"; then

        MPI_FC=${FC}

      else
        MPI_FC_CANDIDATES=""
    
        if test -n "${FC}"; then
          MPI_FC_CANDIDATES=${FC}
        elif test -n "${FTN}"; then
          MPI_FC_CANDIDATES=${FTN}
        elif test -n "${F90}"; then
          MPI_FC_CANDIDATES=${F90}
        fi

        if test -n "${MPI_FC_CANDIDATES}" && test -f ${MPI_FC_CANDIDATES} ; then
          MPI_FC=${MPI_FC_CANDIDATES}
        else
          if test -z "${MPI_FC_CANDIDATES}"; then
            MPI_FC_CANDIDATES="[mpif90 mpif77]"
          fi

          AC_PATH_PROGS(MPI_FC, ${MPI_FC_CANDIDATES}, [notFound], [PATH = ${MPI_SEEK_PATH}])

          if test "${MPI_FC}" != "notFound" ; then
            FC=${MPI_FC}
          else
            echo "-----"
            echo "Cannot find MPI Fortran 90 compiler in " ${MPI_SEEK_PATH}
            echo "Specify a path to all mpi compilers with --with-mpi-compilers=PATH,"
            echo "or specify a path to top mpi directory (above bin) with --with-mpi=PATH,"
            echo "or specify a fortran 90 compiler using FC=<compiler>"
            echo "-----"
            AC_MSG_ERROR([MPI C compiler not found.])
          fi
        fi
      fi
    fi
  fi
fi

#AC_ARG_WITH(mpi-include,
#[AC_HELP_STRING([--with-mpi-include],[Obsolete.  Use --with-mpi-incdir=DIR instead.  Do not prefix DIR with '-I'.])],
#[AC_MSG_ERROR([--with-mpi-include is an obsolte option.   Use --with-mpi-incdir=DIR instead.  Do not prefix DIR with '-I'.  For example '--with-mpi-incdir=/usr/lam_path/include'.])]
#)

AC_ARG_WITH(mpi-libs,
[AC_HELP_STRING([--with-mpi-libs="LIBS"],[MPI libraries @<:@"-lmpi"@:>@])],
[
  MPI_LIBS=${withval}
  AC_MSG_CHECKING(user-defined MPI libraries)
  AC_MSG_RESULT([${MPI_LIBS}])
]
)

AC_ARG_WITH(mpi-incdir,
[AC_HELP_STRING([--with-mpi-incdir=DIR],[MPI include directory @<:@MPIROOT/include@:>@  Do not use -I])],
[
  MPI_INC=${withval}
  AC_MSG_CHECKING(user-defined MPI includes)
  AC_MSG_RESULT([${MPI_INC}])
]
)

AC_ARG_WITH(mpi-libdir,
[AC_HELP_STRING([--with-mpi-libdir=DIR],[MPI library directory @<:@MPIROOT/lib@:>@  Do not use -L])],
[
  MPI_LIBDIR=${withval}
  AC_MSG_CHECKING(user-defined MPI library directory)
  AC_MSG_RESULT([${MPI_LIBDIR}])
]
)

AC_MSG_CHECKING(whether we are using MPI)
AC_MSG_RESULT([${HAVE_PKG_MPI}])

if test "X${HAVE_PKG_MPI}" = "Xyes"; then
   AC_DEFINE(HAVE_MPI,,[define if we want to use MPI])
fi

dnl Define Automake version of HAVE_MPI if appropriate

AM_CONDITIONAL(HAVE_MPI, [test "X${HAVE_PKG_MPI}" = "Xyes"])

])
