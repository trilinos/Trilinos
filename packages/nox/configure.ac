# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
# 
#            NOX: An Object-Oriented Nonlinear Solver Package
#                 Copyright (2002) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#   
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#   
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
# 
# Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
# (rppawlo@sandia.gov).
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# Initialization 
# ------------------------------------------------------------------------

# This must be the first line in configure.ac.
# Optional 3rd argument is email address for bugs.

AC_INIT(nox, 0.1b)

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.

AC_CONFIG_SRCDIR(src/NOX_Solver_Newton.C)

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

AC_CONFIG_AUX_DIR(config)

# Configure should create src/NOX_Config.h from src/NOX_Config.h.in

AM_CONFIG_HEADER(src/NOX_Config.h:src/NOX_Config.h.in)

# Define $build, $host, $target, etc

AC_CANONICAL_TARGET

# Use automake, getting name and version from variables defined by
# AC_INIT command.

AM_INIT_AUTOMAKE(${PACKAGE_NAME}, ${PACKAGE_VERSION})

# Specify required version of autoconf.

AC_PREREQ(2.52)

# Specify revision number of configure.ac that was used to generate
# configure.

dnl Commented out this line for now because it causes headaches in
dnl committing.
dnl AC_REVISION($Revision$)

# ------------------------------------------------------------------------
# Check for user overrides (to be used later)
# ------------------------------------------------------------------------

if test -n "${CXXFLAGS}"; then
  USERDEF_CXXFLAGS=yes
fi

# ------------------------------------------------------------------------
# MPI configure command line option
# ------------------------------------------------------------------------

AC_ARG_ENABLE(mpi,
[  --enable-mpi            enable MPI using mpiCC C++ compiler],
[
HAVE_PKG_MPI=yes 
MPI_CXX=mpiCC
],
[HAVE_PKG_MPI=no]
)

AC_ARG_WITH(mpi,
[  --with-mpi              specify root directory of MPI installation],
[
HAVE_PKG_MPI=yes
MPI_DIR=${withval}
AC_MSG_CHECKING(MPI directory)
AC_MSG_RESULT([${MPI_DIR}])
]
)

AC_ARG_WITH(mpi-cxx,
[  --with-mpi-cxx          specify MPI C++ compiler],
[
HAVE_PKG_MPI=yes 
MPI_CXX=${withval}
AC_MSG_CHECKING(user-defined MPI C++ compiler)
AC_MSG_RESULT([${MPI_CXX}])
]
)

AC_ARG_WITH(mpi-include,
[  --with-mpi-include      specify include directory for MPI],
[
HAVE_PKG_MPI=yes
MPI_INC=${withval}
AC_MSG_CHECKING(user-defined MPI includes)
AC_MSG_RESULT([${MPI_INC}])
]
)

AC_ARG_WITH(mpi-libs,
[  --with-mpi-libs         specify MPI libraries],
[
HAVE_PKG_MPI=yes
MPI_LIBS=${withval}
AC_MSG_CHECKING(user-defined MPI libraries)
AC_MSG_RESULT([${MPI_LIBS}])
]
)

AC_ARG_WITH(mpi-libdir,
[  --with-mpi-libdir       specify location of MPI libraries],
[
HAVE_PKG_MPI=yes
MPI_LIBDIR=${withval}
AC_MSG_CHECKING(user-defined MPI libraries)
AC_MSG_RESULT([${MPI_LIBS}])
]
)

AC_MSG_CHECKING(whether we are using MPI)
AC_MSG_RESULT([${HAVE_PKG_MPI}])

# ------------------------------------------------------------------------
# Check for MPI compiler (must be done *before* AC_PROG_CXX)
# ------------------------------------------------------------------------

if test -n "${MPI_CXX}"; then
  AC_CHECK_PROG(MPI_CXX_EXISTS, ${MPI_CXX}, yes, no)
  if test "X${MPI_CXX_EXISTS}" = "Xyes"; then
    CXX=${MPI_CXX}
  else
    AC_MSG_ERROR([MPI C++ compiler (${MPI_CXX}) not found.])
  fi
fi

# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

AC_PROG_CXX(CC c++ g++ gcc cxx)
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Checks for extra flags, etc
# ------------------------------------------------------------------------

AC_ARG_WITH(flags,
[  --with-flags            add extra preprocessor flags ("-Dflag -Idir")],
[XTRAFLAGS=${withval}]
)

AC_MSG_CHECKING([for extra preprocessor flags])
if test -z "${XTRAFLAGS}"; then
  AC_MSG_RESULT(none)
else
  AC_MSG_RESULT([${XTRAFLAGS}])
  CPPFLAGS="${CPPFLAGS} ${XTRAFLAGS}"
fi

# ------------------------------------------------------------------------
# Checks for special compiler options
# ------------------------------------------------------------------------

AC_ARG_ENABLE(debug,
[  --enable-debug          compile with debugging],
[DEBUG=${enableval}],
[DEBUG="default"]
)
AC_MSG_CHECKING(whether debugging is on)
AC_MSG_RESULT([${DEBUG}])

AC_ARG_ENABLE(opt,
[  --enable-opt            compile with optimization],
[OPT=${enableval}],
[OPT=default]
)
AC_MSG_CHECKING([whether optimization is on])
AC_MSG_RESULT([${OPT}])

case $target in

  # IBM test
  rs6000-ibm-aix3.2.* | rs6000-ibm-aix4.* | powerpc-ibm-aix4.*)

    case $CXX in
      xlC* ) SPECIAL_CXXFLAGS="-qnolm -qrtti";;
    esac

  ;;

  # SGI/IRIX test
  mips-sgi-irix* ) 

    if test "${CXX}" = "CC"; then
     SPECIAL_CXXFLAGS="-LANG:std -LANG:ansi-for-init-scope=ON -64 -DSGI64"
     SPECIAL_AR="CC -ar -W,l -W,s -o"    
    fi

  ;;

esac

AC_MSG_CHECKING([for special C++ flags])
if test "X${USERDEF_CXXFLAGS}" = "Xyes"; then

  AC_MSG_RESULT([${CXXFLAGS} (user-specified via CXXFLAGS)])

else

  if test -n "${SPECIAL_CXXFLAGS}"; then
    AC_MSG_RESULT([${SPECIAL_CXXFLAGS}])
    CXXFLAGS="${CXXFLAGS} ${CXXOPTIONS}"
  else
    AC_MSG_RESULT([${CXXFLAGS} (default)])
  fi

fi

# ------------------------------------------------------------------------
# Alternate archiver
# ------------------------------------------------------------------------

AC_ARG_WITH(ar,
[  --with-ar               override archiver command (default is "ar cru")],
[
AC_MSG_CHECKING(user-defined archiver)
AC_MSG_RESULT([${withval}])
USE_ALTERNATE_AR=yes
ALTERNATE_AR="${withval}"
]
)

if test -n "${SPECIAL_AR}" && test "X${USE_ALTERNATE_AR}" != "Xyes"; then
  USE_ALTERNATE_AR=yes
  ALTERNATE_AR="${SPECIAL_AR}"
fi

AC_MSG_CHECKING(for special archiver command)
if test "X${USE_ALTERNATE_AR}" = "Xyes"; then
   AC_MSG_RESULT([${ALTERNATE_AR}])
   AM_CONDITIONAL(USE_ALTERNATE_AR, true)
else
   AC_MSG_RESULT([none])
   AM_CONDITIONAL(USE_ALTERNATE_AR, false)
fi
AC_SUBST(ALTERNATE_AR)



# ------------------------------------------------------------------------
# MPI link check
# ------------------------------------------------------------------------

if test "X${HAVE_PKG_MPI}" = "Xyes"; then

  if test -n "${MPI_DIR}" && test -z "${MPI_INC}"; then
    MPI_INC=${MPI_DIR}/include
  fi

  if test -n "${MPI_DIR}" && test -z "${MPI_LIBDIR}"; then
    MPI_LIBDIR=${MPI_DIR}/lib
  fi

  if test -n "${MPI_INC}"; then
    CPPFLAGS="${CPPFLAGS} ${MPI_INC}"
  fi

  if test -n "${MPI_LIBDIR}"; then
    LDFLAGS="${LDFLAGS} ${MPI_LIBDIR}"
  fi

  if test -z "${MPI_LIBS}" && test -n "${MPI_LIBDIR}"; then
    MPI_LIBS="-lmpi"
  fi

  if test -n "${MPI_LIBS}"; then
    LIBS = "${MPI_LIBS} ${LIBS}"
  fi

  AC_LANG_CPLUSPLUS 
  AC_MSG_CHECKING(whether MPI will link using C++ compiler)
  AC_TRY_LINK([#include <mpi.h>],
  [int c; char** v; MPI_Init(&c,&v);],
  [AC_MSG_RESULT(yes)], 
  [AC_MSG_RESULT(no)  
   AC_MSG_ERROR(MPI cannot link)]
  )

fi

# ------------------------------------------------------------------------
# Checks for header files
# ------------------------------------------------------------------------

#iostream
AC_LANG_CPLUSPLUS 
AC_MSG_CHECKING([for iostream])
AC_TRY_LINK([#include <iostream>],[],
[
  AC_MSG_RESULT(yes)
  AC_DEFINE(HAVE_IOSTREAM,,[define if we can include <iostream>])
], 
[
  AC_MSG_RESULT(no)
  AC_MSG_CHECKING([for iostream.h])
  AC_TRY_LINK([#include <iostream.h>],[],
  [
    AC_MSG_RESULT(yes)
    AC_DEFINE(HAVE_IOSTREAM_H,,[define if we can include <iostream.h>])
  ], 
  [
    AC_MSG_RESULT(no)
    AC_MSG_NOTICE([You do not appear to have iostream nor iostream.h!!])
  ]
  )
]
)

#iomanip
AC_LANG_CPLUSPLUS 
AC_MSG_CHECKING([for iomanip])
AC_TRY_LINK([#include <iomanip>],[],
[
  AC_MSG_RESULT(yes)
  AC_DEFINE(HAVE_IOMANIP,,[define if we can include <iomanip>])
], 
[
  AC_MSG_RESULT(no)
  AC_MSG_CHECKING([for iomanip.h])
  AC_TRY_LINK([#include <iomanip.h>],[],
  [
    AC_MSG_RESULT(yes)
    AC_DEFINE(HAVE_IOMANIP_H,,[define if we can include <iomanip.h>])
  ], 
  [
    AC_MSG_RESULT(no)
    AC_MSG_NOTICE([You do not appear to have iomanip nor iomanip.h!!])
  ]
  )
]
)

#cmath
AC_LANG_CPLUSPLUS 
AC_MSG_CHECKING([for cmath])
AC_TRY_LINK([#include <cmath>],[],
[
  AC_MSG_RESULT(yes)
  AC_DEFINE(HAVE_CMATH,,[define if we can include <cmath>])
], 
[
  AC_MSG_RESULT(no)
  AC_MSG_CHECKING([for math.h])
  AC_TRY_LINK([#include <math.h>],[],
  [
    AC_MSG_RESULT(yes)
    AC_DEFINE(HAVE_MATH_H,,[define if we can include <math.h>])
  ], 
  [
    AC_MSG_RESULT(no)
    AC_MSG_NOTICE([You do not appear to have cmath nor math.h!!])
  ]
  )
]
)

# can we say using std::ostream?
AC_LANG_CPLUSPLUS 
AC_MSG_CHECKING(whether "using std::ostream" is okay)
AC_TRY_LINK(
[
#ifdef HAVE_IOSTREAM
#include <iostream>
#elif defined(HAVE_IOSTREAM_H)
#include <iostream.h>
#endif
],
[using std::ostream;],
[AC_MSG_RESULT(yes)], 
[
  AC_MSG_RESULT(no)  
  AC_DEFINE(NO_USING_STD_OSTREAM,,[define if "using std::ostream" is okay])
]
)


# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))

# ------------------------------------------------------------------------
# Specify package directory
# ------------------------------------------------------------------------

AC_ARG_WITH(pkgdir,
[  --with-pkgdir           specify root directory of add'l packages],
[
USE_PKGDIR=yes
PKGDIR=${withval}
]
)

if test "X${USE_PKGDIR}" = "Xyes"; then

   AC_MSG_CHECKING(whether package directory is valid)
   if test -d ${PKGDIR}; then
      AC_MSG_RESULT(yes)
   else
      AC_MSG_NOTICE(no)
      AC_MSG_ERROR([Specified package directory does not exist])
   fi

   LDFLAGS="${LDFLAGS} -L${PKGDIR}/lib"
   CPPFLAGS="${CPPFLAGS} -I${PKGDIR}/include"

fi

# ------------------------------------------------------------------------
# NOXEPETRA LIBRARY
# ------------------------------------------------------------------------

AC_ARG_ENABLE(noxepetra,
[  --disable-noxepetra     disable building noxepetra library],
[BUILD_NOXEPETRA=${enableval}],
[BUILD_NOXEPETRA=yes]
)

AC_ARG_ENABLE(trilinos-dirs,
[  --enable-trilinos-dirs  enable use of Trilinos directory structure],
[USE_TRILINOS_DIRS=${withval}],
[USE_TRILINOS_DIRS=default]
)

AC_DEFUN(TRILINOS_LIB_CHECK,
[
  AC_ARG_WITH($1-libdir,
  [  --with-$1-libdir    specify location of $1 library],
  [$1_LIBDIR=$withval])

  AC_MSG_CHECKING(for $1)
  if test -n "${$1_LIBDIR}" && test -f ${$1_LIBDIR}/lib$1.a; then
    LDFLAGS="-L${$1_LIBDIR} ${LDFLAGS}"
    AC_MSG_RESULT(yes)
  elif test -f ${TRILINOS_LIBDIR}/lib$1.a; then
    AC_MSG_RESULT(yes)
  elif test -f ${TRILINOS_LIBDIR}/$1/lib$1.a; then
    LDFLAGS="-L${TRILINOS_LIBDIR}/$1 ${LDFLAGS}"
    AC_MSG_RESULT(yes)
  else
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([Cannot find $1 library])
  fi  

  LIBS="-l$1 ${LIBS}"

])


# Case I: Using standard Trilinos directory structure
# (a) Check if TRILINOS_HOME is valid
if test "X${BUILD_NOXEPETRA}" = "Xyes"; then

  dnl Checks for libraries.
  AC_F77_LIBRARY_LDFLAGS
  LIBS="${LIBS} ${FLIBS}"

  # Find valid BLAS
  AC_LANG_FORTRAN77
  AC_SEARCH_LIBS(dcopy,[cxml blas],,AC_MSG_ERROR([Cannot find BLAS library]))

  # Find valid LAPACK
  AC_LANG_FORTRAN77
  AC_SEARCH_LIBS(dgetrs,[cxml lapack],,AC_MSG_ERROR([Cannot find LAPACK library]))

  # Set TRILINOS_COMM (used several places)
  if test "X${HAVE_PKG_MPI}" = "Xyes"; then
    TRILINOS_COMM="MPI"
  else
    TRILINOS_COMM="SERIAL"
  fi

  # If the package directory isn't specified, use the standard
  # TRILINOS directory structure
  if test "X${USE_TRILINOS_DIRS}" != "Xno" && test "X${USE_PKGDIR}" != "Xyes"; then
    USE_TRILINOS_DIRS=yes
  fi

  if test "X${USE_TRILINOS_DIRS}" != "Xyes" && test "X${USE_PKGDIR}" != "Xyes"; then
    AC_MSG_ERROR([Must use --enable-trilinos-dirs or --with-pkgdir])
  fi

  if test "X${USE_TRILINOS_DIRS}" = "Xyes"; then

    AC_MSG_CHECKING(whether TRILINOS_HOME is defined)
    if test -n "${TRILINOS_HOME}"; then
       AC_MSG_RESULT([yes (${TRILINOS_HOME})])
    else
       AC_MSG_RESULT([no (using ${HOME}/Trilinos)])
       TRILINOS_HOME=${HOME}/Trilinos
    fi
 
    AC_MSG_CHECKING(whether TRILINOS_HOME is valid)
    if test -d ${TRILINOS_HOME}; then
       AC_MSG_RESULT(yes)
    else
       AC_MSG_RESULT(no)
       AC_MSG_ERROR([TRILINOS_HOME is invalid])
    fi
 
    AC_MSG_CHECKING(whether TRILINOS_ARCH is defined)
    if test -n "${TRILINOS_ARCH}"; then
       AC_MSG_RESULT([yes (${TRILINOS_ARCH})])
    else
       case $target in
 	rs6000-ibm-aix3.2.* | rs6000-ibm-aix4.* | powerpc-ibm-aix4.*)
   	  TRILINOS_ARCH="AIX"
 	  ;;
 	alpha-unknown-linux-gnu)
 	  TRILINOS_ARCH="LINUXALPHA"
 	  ;;	
 	*)
 	  TRILINOS_ARCH="LINUX"
   	  ;;
       esac
       AC_MSG_RESULT([no (using ${TRILINOS_ARCH})])
    fi
 
 
    AC_MSG_CHECKING(whether TRILINOS_ID is defined)
    if test -n "${TRILINOS_ID}"; then
       AC_MSG_RESULT([yes (${TRILINOS_ID})])
    else
       AC_MSG_RESULT(no)
    fi
 
    AC_MSG_CHECKING(whether TRILINOS_TARGET is defined)
    if test -n "${TRILINOS_TARGET}"; then
       AC_MSG_RESULT([yes (${TRILINOS_TARGET})])
    else
       TRILINOS_TARGET="${TRILINOS_ARCH}.${TRILINOS_COMM}${TRILINOS_ID}"
       AC_MSG_RESULT([no (using ${TRILINOS_TARGET})])
    fi

    EPETRA_INC=${TRILINOS_HOME}/packages/epetra/src
    AZTECOO_INC=${TRILINOS_HOME}/packages/aztecoo/src
    TRILINOS_LIBDIR=${TRILINOS_HOME}/lib/${TRILINOS_TARGET}
    LDFLAGS="${LDFLAGS} -L${TRILINOS_LIBDIR}"

  else

    EPETRA_INC=${PKGDIR}/include/epetra
    AZTECOO_INC=${PKGDIR}/include/aztecoo
    TRILINOS_LIBDIR=${PKGDIR}/lib

  fi

  # Find a valid y12m
  AC_LANG_FORTRAN77
  AC_SEARCH_LIBS(y12mae,[y12m],,AC_MSG_ERROR([Cannot find y12m library]))

  # Check for ifpack
  TRILINOS_LIB_CHECK(ifpack)
  TRILINOS_LIB_CHECK(aztecoo)
  TRILINOS_LIB_CHECK(epetra)

  # If we make it through all the tests, define the correct libs and includes!!
  TRILINOS_CXXFLAGS="-DEPETRA_${TRILINOS_COMM} -I${EPETRA_INC} -I${AZTECOO_INC}"
fi

AC_SUBST(TRILINOS_CXXFLAGS)

# ------------------------------------------------------------------------
# Checks for TRILINOS status 
# ------------------------------------------------------------------------

if test "X${BUILD_NOXEPETRA}" = "Xyes"; then
   AM_CONDITIONAL(BUILD_NOXEPETRA, true)
else
   AM_CONDITIONAL(BUILD_NOXEPETRA, false)
fi

# ------------------------------------------------------------------------
# Checks for types
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for structures
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

# THIS TEST ISN'T WORKING - WILL COME BACK TO IT LATER.
# Check that all libraries will link with C++
dnl if test X${LIBS} != X; then
dnl   AC_MSG_CHECKING([[whether libraries (${LIBS}) will link]])
dnl   AC_LANG_CPLUSPLUS 
dnl   AC_TRY_LINK([],[],
dnl   [AC_MSG_RESULT(yes)],
dnl   [
dnl   AC_MSG_RESULT(no)
dnl   AC_MSG_ERROR(cannot link)
dnl   ]
dnl   )
dnl fi

# ------------------------------------------------------------------------
# Checks for compiler characteristics
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for linker characteristics
# ------------------------------------------------------------------------

# -- DEC CXX --
dnl if test X${target_cpu} = Xalpha && test X${CXX} = Xcxx; then
dnl    XTRALDADD="cxx_repository/*.o"
dnl fi

dnl AC_MSG_CHECKING([for extra object files to be linked])
dnl if test X$XTRALDADD = X; then
dnl AC_MSG_RESULT(none)
dnl else
dnl AC_MSG_RESULT($XTRALDADD)
dnl fi

AC_SUBST(XTRALDADD)

# ------------------------------------------------------------------------
# Checks for library functions
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for system services
# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
# Perform substitutions in output files
# ------------------------------------------------------------------------

AC_SUBST(ac_aux_dir)

# ------------------------------------------------------------------------
# Output files
# ------------------------------------------------------------------------

AC_CONFIG_FILES([
		Makefile 
		src/Makefile 
		example/Makefile 
		example/epetra_DS6.5.1/Makefile
		example/epetra_1Dfem_linear/Makefile
		example/epetra_1Dfem_nonlinear/Makefile
		])

AC_OUTPUT()


