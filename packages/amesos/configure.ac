# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------
#
# Refer to configuration for more information about adding a new 
# package to Amesos
#
# @HEADER
# ************************************************************************
# 
#                 Amesos: Direct Sparse Solver Package
#                 Copyright (2004) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
# @HEADER

#  Notes:  Ken Stanley 23 JAN 2004 changed AC_CHECK_FILE to TAC_CHECK_FILE
#  to allow cofigure to work on cross-compiled systems.  This is far from 
#  ideal as many files ought to be checked on cross-compiled systems. 

# ------------------------------------------------------------------------
# Initialization 
# ------------------------------------------------------------------------

# This must be the first line in configure.ac.
# Optional 3rd argument is email address for bugs.

AC_INIT([amesos],[0.1],[ksstanl@sandia.gov])

# Hello World!
echo "----------------------------------------"
echo "Running Amesos Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.

AC_CONFIG_SRCDIR(/src/Amesos_Dscpack.cpp)

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

AC_CONFIG_AUX_DIR(config)

# Configure should create src/Amesos_config.h from src/Amesos_config.h.in

AM_CONFIG_HEADER(src/Amesos_config.h:src/Amesos_config.h.in)

# Allow users to specify their own "install" command.  If none is specified,
# the default is install-sh found in the config subdirectory.
                                                                                
AC_ARG_WITH(install,
 [AC_HELP_STRING([--with-install=INSTALL_PROGRAM],
 [Use the installation program INSTALL_PROGRAM rather the default that is provided.  For example --with-install="/path/install -p"])],
 [
   INSTALL=$withval
   INSTALL_PROGRAM=$withval
   INSTALL_SCRIPT=$withval
   INSTALL_DATA="$withval -m 644"
 ],)
                                                                                
# AM_MAINTAINER_MODE turns off maintainer-only makefile targets by
# default, and changes configure to understand a
# --enable-maintainer-mode option. --enable-maintainer-mode turns the
# maintainer-only targets back on. The maintainer-only makefile
# targets permit end users to clean automatically-generated files such
# as configure, which means they have to have autoconf and automake
# installed to repair the damage. AM_MAINTAINER_MODE makes it a bit
# harder for users to shoot themselves in the foot.

AM_MAINTAINER_MODE

# Define $build, $host, $target, etc

AC_CANONICAL_TARGET

# Use automake

AM_INIT_AUTOMAKE(1.8.2)

# Specify required version of autoconf.

AC_PREREQ(2.59)

# Specify revision number of configure.ac that was used to generate
# configure.

dnl Commented out this line for now because it causes headaches in
dnl committing.
dnl AC_REVISION($Revision$)

# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI


# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ gcc c++ cxx)
BASE_CXX=${CXX}
AC_PROG_F77(f77 g77 f90 xlf90)
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Checks for --enable-debug, --enable-optimize, --enable-profile, --enable-purify
# Changes flags and compiler names accordingly
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_COMPILER_OPTIONS

# ------------------------------------------------------------------------
# Checks for extra flags
# ------------------------------------------------------------------------
# Check if --with-flags present, prepend any specs to FLAGS

TAC_ARG_WITH_FLAGS(ccflags, CCFLAGS)
TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
# No C or Fortran code needs to be compiled at this time
TAC_ARG_WITH_FLAGS(cflags, CFLAGS)
TAC_ARG_WITH_FLAGS(fflags, FFLAGS)
TAC_ARG_WITH_FLAGS(ldflags, LDFLAGS)

# ------------------------------------------------------------------------
# Alternate archiver
# ------------------------------------------------------------------------

TAC_ARG_WITH_AR

# ------------------------------------------------------------------------
# MPI link check
# ------------------------------------------------------------------------
TAC_ARG_CHECK_MPI

# ------------------------------------------------------------------------
# Special check for case of running on TFLOP - not needed for now
# ------------------------------------------------------------------------
#AC_MSG_CHECKING([if host is TFLOP (and therefore needs special TFLOP define statement)])
#if test "X${host}" = "Xtflop"; then
#  AC_MSG_RESULT(yes)
#  AC_DEFINE(TFLOP,,[define if compiling for ASCI Red, aka TFLOP])
#else
#  AC_MSG_RESULT(no)
#fi

# ------------------------------------------------------------------------
# Checks for header files
# ------------------------------------------------------------------------

AC_CHECK_HEADERS([cstdlib stdlib.h cstdio stdio.h cassert assert.h],,,)
AC_CHECK_HEADERS([iomanip iomanip.h iostream iostream.h string string.h],,,)

# **AMESOS SPECIFIC**
# The next two rows of headers are not verified to be present below.  At this time
# configure could complete without required headers.  We need to upgrade to include
# c std and vanilla versions of these headers where appropriate.  A centralized ConfigDefs
# file should also be written.  Errors could also occur if the c std version of a header is
# present because currently the vanilla header is simply included - no check is done yet.

AC_CHECK_HEADERS([vector map time.h fstream sys/utsname.h unistd.h],,,)
AC_CHECK_HEADERS([sys/time.h sys/resource.h],,,)

if test "X$ac_cv_header_cstdlib" != "Xyes" && test "X$ac_cv_header_stdlib_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdlib or stdlib.h.  Amesos requires at least one or the other.  See http://software.sandia.gov/trilinos/packages/amesos/ for details]])
fi

if test "X$ac_cv_header_cstdio" != "Xyes" && test "X$ac_cv_header_stdio_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdio or stdio.h.  Amesos requires at least one or the other.]])
fi

if test "X$ac_cv_header_cassert" != "Xyes" && test "X$ac_cv_header_assert_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cassert or assert.h.  Amesos requires at least one or the other.]])
fi

if test "X$ac_cv_header_string" != "Xyes" && test "X$ac_cv_header_string_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find string or string.h.  Amesos requires at least one or the other.]])
fi

if test "X$ac_cv_header_iostream" != "Xyes" && test "X$ac_cv_header_iostream_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find iostream or iostream.h.  Amesos requires at least one or the other.]])
fi

if test "X$ac_cv_header_iomanip" != "Xyes" && test "X$ac_cv_header_iomanip_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find iomanip or iomanip.h.  Amesos requires at least one or the other.]])
fi  

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE

# ------------------------------------------------------------------------
# Checks for compiler characteristics
# ------------------------------------------------------------------------
AC_CXX_NAMESPACES
AC_CXX_HAVE_STL
AC_CXX_BOOL
AC_CXX_MUTABLE
AC_CXX_NEW_FOR_SCOPING
AC_CXX_STD_SPRINTF

# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------

#TAC_ARG_ENABLE_OPTION(threads, [enable shared memory threads], THREADS, no)
#AM_CONDITIONAL(HAVE_THREADS, [test "X$ac_cv_use_threads" != "Xno"])

#TAC_ARG_ENABLE_OPTION(fatal_messages, [enable printing of fatal error messages], 
#                   FATAL_MESSAGES, yes)
#TAC_ARG_ENABLE_OPTION(warning_messages, [enable printing of warning error messages], 
#                   WARNING_MESSAGES, no)
#TAC_ARG_ENABLE_OPTION(format_io, [Use package-defined IO formatting for output], 
#                   FORMAT_IO, yes)
#TAC_ARG_WITH_PACKAGE(zoltan, [Enable Zoltan interface support], ZOLTAN, no)
#TAC_ARG_WITH_PACKAGE(rtop, [Enable RTOp interface support], RTOP, no)
#TAC_ARG_WITH_PACKAGE(esi, [Enable ESI interface support], ESI, yes)
#TAC_ARG_WITH_PACKAGE(superlu, [Enable SuperLU interface support], SUPERLU, no)

# ------------------------------------------------------------------------
# Checks if tests and examples should be built
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE(tests, [Build tests for all Trilinos packages (not all packages are sensitive to this option)], TESTS, yes)
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_tests" != "Xno")

TAC_ARG_ENABLE_FEATURE(examples, [Build examples for all Trilinos packages (not all packages are sensitive to this option)], EXAMPLES, yes)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_examples" != "Xno")

# ------------------------------------------------------------------------
# Specify other directories 
# ------------------------------------------------------------------------

# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS

# Define F77_FUNC that will be used to link with Fortran subroutines, e.g., BLAS/LAPACK
AC_F77_WRAPPERS

# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

  # Find valid SuperLU Dist and Kundert libraries
### In the future we need to be able to disable dependence on Kundert lib

## if test "X$ADDON_SLUD" != "Xno"; then
##   AC_SEARCH_LIBS(pdgssvx,[superludist],,
##   AC_MSG_ERROR([[Cannot find SuperLU Dist library, specify a location using --with-ldflags (ex. --with-ldflags=-L/usr/path/lib).  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]]))  
## fi

dnl Replace `main' with a function in -lm:
#AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
#AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))

dnl Check for user-specified names for LAPACK and BLAS
TAC_ARG_WITH_BLASLIB
# ------------------------------  ABC
TAC_ARG_WITH_LAPACKLIB
# ------------------------------  ABCD

  ACX_BLAS([], [AC_MSG_ERROR([[Cannot find blas library, specify a path using --with-ldflags=-I<DIR> (ex. --with-ldflags=-I/usr/path/lib) or a specific library using --with-blas=DIR/LIB (for example --with-blas=/usr/path/lib/libcxml.a).  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]])])
# ------------------------------  ABCDE
  ACX_LAPACK([], [AC_MSG_ERROR([[Cannot find lapack library, specify a path using --with-ldflags=-I<DIR> (ex. --with-dlflags=-I/usr/path/lib) or a specific library using --with-lapack=DIR/LIB (ex. --with-lapack=/usr/path/lib/libcxml.a).  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]])])

# ------------------------------------------------------------------------
# Checks for addons
# ------------------------------------------------------------------------

TAC_ARG_WITH( trilinos3pldir, [Set Trilinos3PL directory], must_specify_trilinos3pldir )

## echo "  tac_with_trilinos3pldir =  ${tac_with_trilinos3pldir} " 
TAC_CHECK_FILE(${tac_with_trilinos3pldir})

#
#  --enable-amesos-slud, --with-amesos-sluddir, --with-amesos-sludlib, 
#    --with-amesos-sludincdir
#
TAC_ARG_ENABLE_FEATURE_SUB( amesos, slud, [Enable Amesos interface support for SuperLU Dist], 
		AMESOS_SLUD, no)

TAC_ARG_ENABLE_FEATURE_SUB( amesos, superludist, [Enable Amesos interface support for SuperLU Dist], 
		AMESOS_SUPERLUDIST, no)

TAC_ARG_ENABLE_FEATURE_SUB( amesos, slud2, [Enable Amesos interface support for SuperLU Dist 2.0 (requires amesos-slud)], 
		AMESOS_SLUD2, no)

TAC_ARG_ENABLE_FEATURE_SUB( amesos, dscpack, [Enable Amesos interface support for DSCPACK], 
		AMESOS_DSCPACK, no)

TAC_ARG_ENABLE_FEATURE_SUB( amesos, scalapack, [Enable Amesos interface support for SCALAPACK], 
		AMESOS_SCALAPACK, no)

TAC_ARG_ENABLE_FEATURE_SUB( amesos, umfpack, [Enable Amesos interface support for UMFPACK], 
		AMESOS_UMFPACK, no)

TAC_ARG_ENABLE_FEATURE_SUB( amesos, superlu, [Enable Amesos interface support for SUPERLU], 
		AMESOS_SUPERLU, no)

TAC_ARG_ENABLE_FEATURE_SUB( amesos, klu, [Enable Amesos interface support for KLU], 
		AMESOS_KLU, yes)

TAC_ARG_ENABLE_FEATURE_SUB( amesos, mumps, [Enable Amesos interface support for MUMPS], 
		AMESOS_MUMPS, no)

TAC_ARG_ENABLE_FEATURE_SUB(amesos, mumps_mpi_c2f, [Enable conversion from C++ MPI commnicator to FORTRAN MPI communicator.], AMESOS_MPI_C2F, yes) 

TAC_ARG_ENABLE_FEATURE_SUB( amesos, slus, [Enable Amesos interface support for SuperLU Dist], 
		AMESOS_SLUS, no)

TAC_ARG_ENABLE_FEATURE_SUB( amesos, kundert, [Enable Amesos interface support for Kundert], 
		AMESOS_KUNDERT, no)

AM_CONDITIONAL(NEED_PARAMETER_LIST, false ) 

#
#  Adding -tweak to the cxx compile flags causes 
#  all template instatiations to be placed in the library
#  and instructs the linker to ignore the multiple definitions
#  that result.  
#  
#  Ken Stanley added this in June 2003 as the lesser of the 
#  known evil ways of dealing with stl template instantiation
#
case "$target" in
  alpha*-dec-osf*)  
    if test $CXX = cxx; then
      CXXFLAGS="-tweak ${CXXFLAGS}"
    fi
  ;;
esac

#
#
#
      CXXFLAGS="  ${CXXFLAGS} -I\${top_srcdir}/amd/src/ -I\${top_srcdir}/src/klu"
      CFLAGS="  ${CFLAGS} -I\${top_srcdir}/amd/src/"


#
#  Test to see if everything is where it ought to be:
#    Make sure that supermatrix.h is in the SuperLU distributed include directory
#    Make sure that pdgssvx is in the SuperLU distributed library
#
amesos_libs=" " 
# AC_LANG_PUSH
#  Atlantis fails unless we do the test links with the C compiler (instead of c++)
AC_LANG([C])

if test "X$ac_cv_use_amesos_superludist" != "Xno"; then

  if test "X${HAVE_PKG_MPI}" != "Xyes"; then
    AC_MSG_ERROR([[Amesos/SuperludistOO requires MPI.  See http://software.sandia.gov/trilinos/packages/amesos/ for details]])
  fi

  TAC_ARG_WITH_SUB( amesos, superludistdir, [Set SuperLU dist directory],  ${tac_with_trilinos3pldir}/SuperLU_DIST_2.0 )
  TAC_CHECK_FILE(${tac_with_amesos_superludistdir})

  TAC_ARG_WITH_SUB( amesos, superludistlib, [Set SuperLU dist library],  ${tac_with_amesos_superludistdir}/libsuperludist.a )
  TAC_CHECK_FILE(${tac_with_amesos_superludistlib})

  TAC_ARG_WITH_SUB( amesos, superludistincdir, [Set SuperLU dist include directory],  ${tac_with_amesos_superludistdir}/SRC )
  TAC_CHECK_FILE(${tac_with_amesos_superludistincdir})

#
# This is a bad way to check for supermatrix.h, it won't pick supermatrix.h up if
# it is in the include path, i.e. /usr/include - see tac_arg_check_mpi for a better way
#


  TAC_CHECK_FILE(${tac_with_amesos_superludistincdir}/supermatrix.h)

  TAC_CHECK_FILE(${tac_with_amesos_superludistincdir}/supermatrix.h,,
    AC_MSG_ERROR([[Cannot find supermatrix.h in ${tac_with_amesos_superludistincdir} , specify a location using --with-amesos-superludistincdir, --with-amesos-superludistdir or --with-trilinos3pldir.  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]]))


##  LDFLAGS="${tac_with_amesos_superludistlib} ${LDFLAGS}"

  CXXFLAGS="-I${tac_with_amesos_superludistincdir} ${CXXFLAGS}"

    SAVE_LIBS="${LIBS}"
    LIBS="${tac_with_amesos_superludistlib} ${BLAS_LIBS} ${MPI_LIBS} ${FLIBS} ${LIBS}"
     AC_TRY_LINK_FUNC(set_default_options,AC_DEFINE(OLD_SUPERLUDIST,,[Define if we are using the old Superludist]) ) 
     AC_TRY_LINK_FUNC(pdgssvx,,
      AC_MSG_ERROR([[Cannot find SuperLU Dist library, specify a location using --with-amesos-superludistlib, (ex. --with-amesos-superludistlib=/home/kstanley/Trilinos3PL/SuperLU_DIST/libsuperludist.a)]]))  
    LIBS="${SAVE_LIBS}"
  amesos_libs="${amesos_libs} ${tac_with_amesos_superludistlib}"

  AM_CONDITIONAL(HAVE_AMESOS_SUPERLUDIST, true ) 
  AM_CONDITIONAL(NEED_PARAMETER_LIST, true ) 
else
  AM_CONDITIONAL(HAVE_AMESOS_SUPERLUDIST, false ) 
fi




#
#  Test to see if everything is where it ought to be:
#    Make sure that supermatrix.h is in the SuperLU distributed include directory
#    Make sure that pdgssvx is in the SuperLU distributed library
#

if test "X$ac_cv_use_amesos_slud" != "Xno"; then

  if test "X${HAVE_PKG_MPI}" != "Xyes"; then
    AC_MSG_ERROR([[Amesos/SuperludistOO requires MPI.  See http://software.sandia.gov/trilinos/packages/amesos/ for details]])
  fi

  TAC_ARG_WITH_SUB( amesos, sluddir, [Set SuperLU dist directory],  ${tac_with_trilinos3pldir}/SuperLU_DIST )
  TAC_CHECK_FILE(${tac_with_amesos_sluddir})

  TAC_ARG_WITH_SUB( amesos, sludlib, [Set SuperLU dist library],  ${tac_with_amesos_sluddir}/libsuperludist.a )
  TAC_CHECK_FILE(${tac_with_amesos_sludlib})

  TAC_ARG_WITH_SUB( amesos, sludincdir, [Set SuperLU dist include directory],  ${tac_with_amesos_sluddir}/SRC )
  TAC_CHECK_FILE(${tac_with_amesos_sludincdir})

#
# This is a bad way to check for supermatrix.h, it won't pick supermatrix.h up if
# it is in the include path, i.e. /usr/include - see tac_arg_check_mpi for a better way
#


  TAC_CHECK_FILE(${tac_with_amesos_sludincdir}/supermatrix.h)

  TAC_CHECK_FILE(${tac_with_amesos_sludincdir}/supermatrix.h,,
    AC_MSG_ERROR([[Cannot find supermatrix.h in ${tac_with_amesos_sludincdir} , specify a location using --with-amesos-sludincdir, --with-amesos-sluddir or --with-trilinos3pldir.  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]]))


##  LDFLAGS="${tac_with_amesos_sludlib} ${LDFLAGS}"

  CXXFLAGS="-I${tac_with_amesos_sludincdir} ${CXXFLAGS}"

    SAVE_LIBS="${LIBS}"
    LIBS="${tac_with_amesos_sludlib} ${BLAS_LIBS} ${MPI_LIBS} ${FLIBS} ${LIBS}"
     AC_TRY_LINK_FUNC(pdgssvx,,
      AC_MSG_ERROR([[Cannot find SuperLU Dist library, specify a location using --with-amesos-sludlib, (ex. --with-amesos-sludlib=/home/kstanley/Trilinos3PL/SuperLU_DIST/libsuperludist.a)]]))  
    LIBS="${SAVE_LIBS}"
  amesos_libs="${amesos_libs} ${tac_with_amesos_sludlib}"

  AM_CONDITIONAL(HAVE_AMESOS_SLUD, true ) 
else
  AM_CONDITIONAL(HAVE_AMESOS_SLUD, false ) 
fi

AC_MSG_CHECKING([ Here we are above the slud2 stuff  ])

if test "X$ac_cv_use_amesos_slud2" != "Xno"; then
  AC_MSG_CHECKING([ Here we are inside the slud2 conditional ])

  if test "X$ac_cv_use_amesos_slud" != "Xno"; then
     AC_MSG_CHECKING([ Here we are inside both the slud2/slud conditional ])
     AM_CONDITIONAL(HAVE_AMESOS_SLUD2, true )
  else 
     AC_MSG_CHECKING([ Here we are inside the slud else clause ])
     AC_MSG_ERROR([Cannot enable slud2 without enabling slud.])
  fi 

else 
    AC_MSG_CHECKING([ Here we are in the slud2 esle clause  ])
   AM_CONDITIONAL(HAVE_AMESOS_SLUD2, false )  
fi
  AC_MSG_CHECKING([ Here we are after the slud2 stuff ])


#
#  --enable-amesos-dscpack, --with-amesos-dscpackdir, --with-amesos-dscpacklib, 
#    --with-amesos-dscpackincdir
#
#
#  Test to see if everything is where it ought to be:
#    Make sure that dscmain.h is in the directory pointed to by amesos-dscpackincdic
#    Make sure that DSC_Clock is in the library pointed to by amesos-dscpacklib
#
if test "X$ac_cv_use_amesos_dscpack" != "Xno"; then

  if test "X${HAVE_PKG_MPI}" != "Xyes"; then
    AC_MSG_ERROR([[Amesos_Dscpack requires MPI.  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]])
  fi

  TAC_ARG_WITH_SUB( amesos, dscpackdir, [Set DSCPACK directory],  ${tac_with_trilinos3pldir}/DSCPACK1.0 )
  TAC_CHECK_FILE(${tac_with_amesos_dscpackdir})

  TAC_ARG_WITH_SUB( amesos, dscpacklib, [Set DSCPACK library],  ${tac_with_amesos_dscpackdir}/DSC_LIB/dsclibdbl.a )
  TAC_CHECK_FILE(${tac_with_amesos_dscpacklib})

  TAC_ARG_WITH_SUB( amesos, dscpackincdir, [Set DSCPACK include directory],  ${tac_with_amesos_dscpackdir}/DSC_LIB )
  TAC_CHECK_FILE(${tac_with_amesos_dscpackincdir})

#
# This is a bad way to check for dscmain.h, it won't pick dscmain.h up if
# it is in the include path, i.e. /usr/include - see tac_arg_check_mpi for a better way
#


  TAC_CHECK_FILE(${tac_with_amesos_dscpackincdir}/dscmain.h)

  TAC_CHECK_FILE(${tac_with_amesos_dscpackincdir}/dscmain.h,,
    AC_MSG_ERROR([[Cannot find dscmain.h in ${tac_with_amesos_dscpackincdir} , specify a location using --with-amesos-dscpackincdir, --with-amesos-dscpackdir or --with-trilinos3pldir.  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]] ))


##  LDFLAGS="${tac_with_amesos_dscpacklib} ${LDFLAGS}"

  CXXFLAGS="-I${tac_with_amesos_dscpackincdir} ${CXXFLAGS}"

    SAVE_LIBS="${LIBS}"
    LIBS="${tac_with_amesos_dscpacklib} ${MPI_LIBS} ${BLAS_LIBS} ${FLIBS} ${LIBS}"
     AC_TRY_LINK_FUNC(DSC_CloseSolver,,
      AC_MSG_ERROR([[Cannot find the DSCPACK library, specify a location using --with-amesos-dscpacklib, (ex. --with-amesos-dscpacklib=/home/kstanley/Trilinos3PL/DSCPACK1.0/dsclibdbl.a).  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]]))  
    LIBS="${SAVE_LIBS}"
  AM_CONDITIONAL(HAVE_AMESOS_DSCPACK, true ) 
  AM_CONDITIONAL(NEED_PARAMETER_LIST, true ) 
  amesos_libs="${amesos_libs} ${tac_with_amesos_dscpacklib}"
else
  AM_CONDITIONAL(HAVE_AMESOS_DSCPACK, false ) 
fi

#
#  --enable-amesos-umfpack, --with-amesos-umfpackdir, --with-amesos-umfpacklib, 
#    --with-amesos-umfpackamdlib, 
#    --with-amesos-umfpackincdir
#
#
#  Test to see if everything is where it ought to be:
#    Make sure that umfmain.h is in the directory pointed to by amesos-umfpackincdic
#    Make sure that amd.h is in the directory pointed to by amesos-umfpackamdincdic
#    Make sure that umfpack_di_default is in the library pointed to by amesos-umfpacklib
#    Make sure that amd_default is in the library pointed to by amesos-umfpackamdlib
#
if test "X$ac_cv_use_amesos_umfpack" != "Xno"; then

  TAC_ARG_WITH_SUB( amesos, umfpackdir, [Set UMFPACK directory],  ${tac_with_trilinos3pldir}/UMFPACKv4.1/ )
  TAC_CHECK_FILE(${tac_with_amesos_umfpackdir})

  TAC_ARG_WITH_SUB( amesos, umfpacklib, [Set UMFPACK library],  ${tac_with_amesos_umfpackdir}UMFPACK/Lib/libumfpack.a )
  TAC_CHECK_FILE(${tac_with_amesos_umfpacklib})

  TAC_ARG_WITH_SUB( amesos, umfpackamdlib, [Set UMFPACK AMD library],  ${tac_with_amesos_umfpackdir}AMD/Lib/libamd.a )
  TAC_CHECK_FILE(${tac_with_amesos_umfpackamdlib})

  TAC_ARG_WITH_SUB( amesos, umfpackincdir, [Set UMFPACK include directory],  ${tac_with_amesos_umfpackdir}UMFPACK/Include )
  TAC_CHECK_FILE(${tac_with_amesos_umfpackincdir})

  TAC_ARG_WITH_SUB( amesos, umfpackamdincdir, [Set UMFPACK AMD include directory],  ${tac_with_amesos_umfpackdir}AMD/Include )
  TAC_CHECK_FILE(${tac_with_amesos_umfpackincdir})

#
# This is a bad way to check for umfpack.h, it won't pick umfpack.h up if
# it is in the include path, i.e. /usr/include - see tac_arg_check_mpi for a better way
#

  TAC_CHECK_FILE(${tac_with_amesos_umfpackincdir}/umfpack.h,,
    AC_MSG_ERROR([[Cannot find umfpack.h in ${tac_with_amesos_umfpack} , specify a location using --with-amesos-umfpackincdir, --with-amesos-umfpackdir or --with-trilinos3pldir.  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]] ))


  TAC_CHECK_FILE(${tac_with_amesos_umfpackamdincdir}/amd.h,,
    AC_MSG_ERROR([[Cannot find umfpack.h in ${tac_with_amesos_umfpack} , specify a location using --with-amesos-umfpackincdir, --with-amesos-umfpackdir or --with-trilinos3pldir.  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]] ))


##  LDFLAGS="${tac_with_amesos_umfpacklib} ${LDFLAGS}"

  CXXFLAGS="-I${tac_with_amesos_umfpackincdir} -I${tac_with_amesos_umfpackamdincdir}   ${CXXFLAGS}"
  CPPFLAGS="-I${tac_with_amesos_umfpackincdir} -I${tac_with_amesos_umfpackamdincdir}   ${CPPFLAGS}"

    SAVE_LIBS="${LIBS}"
    LIBS="${tac_with_amesos_umfpacklib} ${BLAS_LIBS} ${MPI_LIBS} ${FLIBS} ${LIBS}"
     AC_TRY_LINK_FUNC(umfdi_local_search,,
      AC_MSG_ERROR([[Cannot find the UMFPACK library, specify a location using --with-amesos-umfpacklib, (ex. --with-amesos-umfpacklib=/home/kstanley/Trilinos3PL/UMFPACKv4.1/libumfpack.a).  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]]))  

    LIBS="${tac_with_amesos_umfpackamdlib} ${BLAS_LIBS} ${MPI_LIBS} ${FLIBS} ${LIBS}"
     AC_TRY_LINK_FUNC(amd_defaults,,
      AC_MSG_ERROR([[Cannot find the AMD library needed by Umfpack, specify a location using --with-amesos-umfpackamdlib, (ex. --with-amesos-umfpacklib=/home/kstanley/Trilinos3PL/UMFPACKv4.1/AMD/Lib/libamd.a).  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]]))  
    LIBS="${SAVE_LIBS}"
  AM_CONDITIONAL(HAVE_AMESOS_UMFPACK, true ) 
  AM_CONDITIONAL(NEED_PARAMETER_LIST, true ) 
  amesos_libs="${amesos_libs} ${tac_with_amesos_umfpacklib} ${tac_with_amesos_umfpackamdlib}"
else
  AM_CONDITIONAL(HAVE_AMESOS_UMFPACK, false ) 
fi










#
#  --enable-amesos-klu
#
if test "X$ac_cv_use_amesos_klu" != "Xno"; then
  AM_CONDITIONAL(HAVE_AMESOS_KLU, true ) 
  AM_CONDITIONAL(NEED_PARAMETER_LIST, true ) 
  amesos_libs="${amesos_libs} \${top_builddir}/amd/src/libamd.a"
else
  AM_CONDITIONAL(HAVE_AMESOS_KLU, false ) 
fi



#
#  Look for scalapack and blacs libraries 
#    if --enable-amesos-scalapack or --enable-amesos-mumps
#
if test "X$ac_cv_use_amesos_scalapack" != "Xno" || 
   test "X$ac_cv_use_amesos_mumps" != "Xno"; then
      TAC_ARG_WITH_SUB( amesos, scalapackdir, [Set scalapack directory],  ${tac_with_trilinos3pldir}/scalapack )
      TAC_CHECK_FILE(${tac_with_amesos_scalapackdir})

      TAC_ARG_WITH_SUB( amesos, scalapacklib, [Set scalapack library],  ${tac_with_amesos_scalapackdir}/libscalapack.a )
      TAC_CHECK_FILE(${tac_with_amesos_scalapacklib})

      TAC_ARG_WITH_SUB( amesos, blacsdir, [Set blacs directory],  ${tac_with_trilinos3pldir}/BLACS/LIB/ )
      TAC_CHECK_FILE(${tac_with_amesos_blacsdir})

      TAC_ARG_WITH_SUB( amesos, blacslib, [Set blacs library],  ${tac_with_amesos_blacsdir}/libblacs.a )
      TAC_CHECK_FILE(${tac_with_amesos_blacslib})

      TAC_ARG_WITH_SUB( amesos, blacsclib, [Set blacs C library],  ${tac_with_amesos_blacsdir}/libcblacs.a )
      TAC_CHECK_FILE(${tac_with_amesos_blacsclib})

      AC_F77_FUNC(dgebs2d)

    SAVE_LIBS="${LIBS}"
    LIBS="${tac_with_amesos_blacslib}  ${tac_with_amesos_cblacslib} ${BLAS_LIBS} ${MPI_LIBS} ${FLIBS} ${LIBS}"
     AC_TRY_LINK_FUNC($dgebs2d,,
      AC_MSG_ERROR([[Cannot find blacs library, specify a location using --with-amesos-blacslib, (ex. --with-amesos-blacsdir=/home/kstanley/Trilinos3PL/BLACS/LIB).  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]]))  
    LIBS="${SAVE_LIBS}"

      AM_CONDITIONAL(HAVE_AMESOS_SCALAPACK, true ) 
      AM_CONDITIONAL(NEED_PARAMETER_LIST, true ) 
      amesos_libs="${amesos_libs} ${tac_with_amesos_scalapacklib} ${tac_with_amesos_blacslib} ${tac_with_amesos_blacsclib}  ${tac_with_amesos_blacslib} "

      AC_F77_FUNC(pdgetrf)

    SAVE_LIBS="${LIBS}"
    LIBS="${tac_with_amesos_scalapacklib} ${tac_with_amesos_blacslib}  ${tac_with_amesos_cblacslib} ${tac_with_amesos_blacslib} ${BLAS_LIBS} ${MPI_LIBS} ${FLIBS} ${LIBS}"
     AC_TRY_LINK_FUNC($pdgetrf,,
      AC_MSG_ERROR([[Cannot find scalapack library, specify a location using --with-amesos-scalapacklib, (ex. --with-amesos-scalapacklib=/home/kstanley/Trilinos3PL/SuperLU_DIST_2.0/libsuperludist.a).  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]]))  
    LIBS="${SAVE_LIBS}"

else
  AM_CONDITIONAL(HAVE_AMESOS_SCALAPACK, false ) 
fi





#
#  --enable-amesos-mumps, --with-amesos-mumpsdir, --with-amesos-mumpslib, 
#    --with-amesos-mumpspordlib, 
#    --with-amesos-mumpsincdir
#
#
#  Test to see if everything is where it ought to be:
#    Make sure that umfmain.h is in the directory pointed to by amesos-mumpsincdic
#    Make sure that umfpack_di_default is in the library pointed to by amesos-mumpslib
#    Make sure that amd_default is in the library pointed to by amesos-mumpspordlib
#
if test "X$ac_cv_use_amesos_mumps" != "Xno"; then

  TAC_ARG_WITH_SUB( amesos, mumpsdir, [Set MUMPS directory],  ${tac_with_trilinos3pldir}/MUMPS_4.2_beta/ )
  TAC_CHECK_FILE(${tac_with_amesos_mumpsdir})

    TAC_ARG_WITH_SUB( amesos, mumpslib, [Set MUMPS library],  ${tac_with_amesos_mumpsdir}/lib/libdmumps.a )
    TAC_CHECK_FILE(${tac_with_amesos_mumpslib})

    TAC_ARG_WITH_SUB( amesos, mumpspordlib, [Set MUMPS PORD library],  ${tac_with_amesos_mumpsdir}/lib/libpord.a )
    TAC_CHECK_FILE(${tac_with_amesos_mumpspordlib})


  TAC_ARG_WITH_SUB( amesos, mumpsincdir, [Set MUMPS include directory],  ${tac_with_amesos_mumpsdir}include )
  TAC_CHECK_FILE(${tac_with_amesos_mumpsincdir})

#
# This is a bad way to check for dmumps_c.h, it won't pick dmumps_c.h up if
# it is in the include path, i.e. /usr/include - see tac_arg_check_mpi for a better way
#

  TAC_CHECK_FILE(${tac_with_amesos_mumpsincdir}/dmumps_c.h,,
    AC_MSG_ERROR([[Cannot find dmumps_c.h in ${tac_with_amesos_mumps} , specify a location using --with-amesos-mumpsincdir, --with-amesos-mumpsdir or --with-trilinos3pldir.  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]] ))


##  LDFLAGS="${tac_with_amesos_mumpslib} ${LDFLAGS}"

  CXXFLAGS="-I${tac_with_amesos_mumpsincdir} ${CXXFLAGS}"


#
#    The following attempt to confirm that libmumps.a and lidpord.a are the right libraries
#    fails - presumably because we need to link other libraries in (either fortran 
#    scalapack or blacs libraries.)
#

#    SAVE_LIBS="${LIBS}"
#      LIBS="${tac_with_amesos_mumpslib} ${MPI_LIBS} ${LIBS}"
#      AC_TRY_LINK_FUNC(dmumps_164_,,
#      AC_MSG_ERROR([[Cannot find the MUMPS library, specify a location using --with-amesos-mumpslib, (ex. --with-amesos-mumpslib=/home/kstanley/Trilinos3PL/MUMPS_4.2_beta/lib/libdmumps.a).  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]]))  
#    LIBS="${tac_with_amesos_mumpspordlib} ${MPI_LIBS} ${LIBS}"
#      AC_TRY_LINK_FUNC(qsortUpInts,,
#      AC_MSG_ERROR([[Cannot find the PORD library needed by MUMPS, specify a location using --with-amesos-mumpspordlib, (ex. --with-amesos-mumpslib=/home/kstanley/Trilinos3PL/MUMPS_4.2_beta/lib/libpord.a).  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]]))  
#    LIBS="${SAVE_LIBS}"
#
#  Notes for running with the Intel fortran 90 compiler (ifort):
#
#    We may need the following to run in serial:
#      ${tac_with_amesos_mumpsmpiseqlib} 
#    amesos_libs="${amesos_libs} ${tac_with_amesos_mumpslib} ${tac_with_amesos_mumpspordlib}  -llapack -lblas -lg2c  -L/opt/intel_fc_80/lib -lifport -lifcoremt  -lunwind"
#
  AM_CONDITIONAL(HAVE_AMESOS_MUMPS, true ) 
  AM_CONDITIONAL(NEED_PARAMETER_LIST, true ) 
  amesos_libs="${amesos_libs} ${tac_with_amesos_mumpslib} ${tac_with_amesos_mumpspordlib} ${tac_with_amesos_mumpslib} ${tac_with_amesos_scalapacklib}"
else
  AM_CONDITIONAL(HAVE_AMESOS_MUMPS, false ) 
fi




#
#  --enable-amesos-slus, --with-amesos-slusdir, --with-amesos-sluslib, 
#    --with-amesos-slusincdir
#
#
#  Test to see if everything is where it ought to be:
#    Make sure that supermatrix.h is in the SuperLU serial include directory
#    Make sure that SuperLU_timer is in the SuperLU serial library
#

if test "X$ac_cv_use_amesos_slus" != "Xno"; then

  TAC_ARG_WITH_SUB( amesos, slusdir, [Set SuperLU serial directory],  ${tac_with_trilinos3pldir}/SuperLU )
  TAC_CHECK_FILE(${tac_with_amesos_slusdir})

  TAC_ARG_WITH_SUB( amesos, sluslib, [Set SuperLU serial library],  ${tac_with_amesos_slusdir}/libsuperlu.a )
  TAC_CHECK_FILE(${tac_with_amesos_sluslib})

  TAC_ARG_WITH_SUB( amesos, slusincdir, [Set SuperLU serial include directory],  ${tac_with_amesos_slusdir}/SRC )
  TAC_CHECK_FILE(${tac_with_amesos_slusincdir})

#
# This is a bad way to check for supermatrix.h, it won't pick supermatrix.h up if
# it is in the include path, i.e. /usr/include - see tac_arg_check_mpi for a better way
#


  TAC_CHECK_FILE(${tac_with_amesos_slusincdir}/supermatrix.h)

  TAC_CHECK_FILE(${tac_with_amesos_slusincdir}/supermatrix.h,,
    AC_MSG_ERROR([[Cannot find supermatrix.h in ${tac_with_amesos_slusincdir} , specify a location using --with-amesos-slusincdir, --with-amesos-slusdir or --with-trilinos3pldir.  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]] ))


##  LDFLAGS="${tac_with_amesos_sluslib} ${LDFLAGS}"

  CXXFLAGS="-I${tac_with_amesos_slusincdir} ${CXXFLAGS}"

    SAVE_LIBS="${LIBS}"
    LIBS="${tac_with_amesos_sluslib} ${BLAS_LIBS} ${MPI_LIBS} ${FLIBS} ${LIBS}"
     AC_TRY_LINK_FUNC(SuperLU_timer_,,
      AC_MSG_ERROR([[Cannot find SuperLU serial library, specify a location using --with-amesos-sluslib, (ex. --with-amesos-sluslib=/home/kstanley/Trilinos3PL/SuperLU/libsuperlu.a).  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]]))  
    LIBS="${SAVE_LIBS}"
  amesos_libs=" ${tac_with_amesos_sluslib} ${amesos_libs}"

  AM_CONDITIONAL(HAVE_AMESOS_SLUS, true ) 
else
  AM_CONDITIONAL(HAVE_AMESOS_SLUS, false ) 
fi




#
#  --enable-amesos-superlu, --with-amesos-superludir, --with-amesos-superlulib, 
#    --with-amesos-superluincdir
#
#
#  Test to see if everything is where it ought to be:
#    Make sure that supermatrix.h is in the SuperLU serial include directory
#    Make sure that SuperLU_timer is in the SuperLU serial library
#

if test "X$ac_cv_use_amesos_superlu" != "Xno"; then

  TAC_ARG_WITH_SUB( amesos, superludir, [Set SuperLU serial directory],  ${tac_with_trilinos3pldir}/SuperLU_3.0 )
  TAC_CHECK_FILE(${tac_with_amesos_superludir})

  TAC_ARG_WITH_SUB( amesos, superlulib, [Set SuperLU serial library],  ${tac_with_amesos_superludir}/libsuperlu.a )
  TAC_CHECK_FILE(${tac_with_amesos_superlulib})

  TAC_ARG_WITH_SUB( amesos, superluincdir, [Set SuperLU serial include directory],  ${tac_with_amesos_superludir}/SRC )
  TAC_CHECK_FILE(${tac_with_amesos_superluincdir})

#
# This is a bad way to check for supermatrix.h, it won't pick supermatrix.h up if
# it is in the include path, i.e. /usr/include - see tac_arg_check_mpi for a better way
#


  TAC_CHECK_FILE(${tac_with_amesos_superluincdir}/supermatrix.h)

  TAC_CHECK_FILE(${tac_with_amesos_superluincdir}/supermatrix.h,,
    AC_MSG_ERROR([[Cannot find supermatrix.h in ${tac_with_amesos_superluincdir} , specify a location using --with-amesos-superluincdir, --with-amesos-superludir or --with-trilinos3pldir.  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]] ))


##  LDFLAGS="${tac_with_amesos_superlulib} ${LDFLAGS}"

  CXXFLAGS="-I${tac_with_amesos_superluincdir} ${CXXFLAGS}"

    SAVE_LIBS="${LIBS}"
    LIBS="${tac_with_amesos_superlulib} ${BLAS_LIBS} ${MPI_LIBS} ${FLIBS} ${LIBS}"
     AC_TRY_LINK_FUNC(SuperLU_timer_,,
      AC_MSG_ERROR([[Cannot find SuperLU serial library, specify a location using --with-amesos-superlulib, (ex. --with-amesos-superlulib=/home/kstanley/Trilinos3PL/SuperLU/libsuperlu.a).  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]]))  
    LIBS="${SAVE_LIBS}"
  amesos_libs=" ${tac_with_amesos_superlulib} ${amesos_libs}"

  AM_CONDITIONAL(HAVE_AMESOS_SUPERLU, true ) 
else
  AM_CONDITIONAL(HAVE_AMESOS_SUPERLU, false ) 
fi


#
#  --enable-amesos-kundert
#
#  NOTE:  Kundert is supported, but not tested as of May 13, 2003 
#  HAVE_AMESOS_KUNDERT is undef'd in Amesos_ConfigDefs.h 
#

if test "X$ac_cv_use_amesos_kundert" != "Xno"; then
  AM_CONDITIONAL(HAVE_AMESOS_KUNDERT, true ) 

  # ------------------------------------------------------------------------
  # Setup Platform for Kundert
  # ------------------------------------------------------------------------

 case "$target" in
  *freebsd*) PLATFORM=freebsd
           ;;
  *linux*) PLATFORM=linux
         ;;
  *sgi*) PLATFORM=sgi
       ;;
  *powerpc-apple*) PLATFORM=ppc
      ;;
  *powerpc-ibm-aix*) PLATFORM=aix
      ;;
  *solaris*) PLATFORM=solaris
      ;;
  i386*mingw32msvc) PLATFORM=mingw
       ;;
  alpha*-dec-osf*) PLATFORM=DEC
       ;;
  *)   if test "x$PLATFORM" = "x"; then
         AC_MSG_NOTICE([Unrecognized target, using generic linux for PLATFORM])
         Xyce_ARCH=linux
         PLATFORM=linux
       fi
       ;;
 esac

 AC_MSG_NOTICE([Adding Platform Define $PLATFORM to CCFLAGS])
 CPPFLAGS="$CPPFLAGS -D$PLATFORM"
else
  AM_CONDITIONAL(HAVE_AMESOS_KUNDERT, false ) 
fi


# AC_LANG_POP

AC_SUBST(AMESOS_LIBS,${amesos_libs}) 


## AC_ARG_ENABLE(amesos-slud,
## [AC_HELP_STRING([--enable-amesos-slud],[compile SuperLU Dist interface])],
## [ADDON_SLUD=$enableval],
## [ADDON_SLUD=yes]
## )
## AM_CONDITIONAL(BUILD_SLUDIST, test X${ADDON_SLUD} = "Xyes")

## AC_ARG_ENABLE(amesos-slus,
## [AC_HELP_STRING([--enable-amesos-slus],[compile SuperLU Serial interface])],
## [ADDON_SLUS=$enableval],
## [ADDON_SLUS=yes]
## )
## AM_CONDITIONAL(BUILD_SLUSERIAL, test X${ADDON_SLUS} = "Xyes")

## AC_ARG_ENABLE(amesos-kundert,
## [AC_HELP_STRING([--enable-amesos-kundert],[compile Kundert Sparse])],
## [ADDON_KS=$enableval],
## [ADDON_KS=no]
## )
## AM_CONDITIONAL(HAVE_AMESOS_KUNDERT, test X${ADDON_KS} = "Xyes")

# ------------------------------------------------------------------------
# Checks for compiler characteristics
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for linker characteristics
# ------------------------------------------------------------------------

# Determine libraries needed for linking with Fortran
AC_F77_LIBRARY_LDFLAGS

# -- DEC CXX --
# MS # Add the following 21-Apr-04
case $target in
  alpha*)
    if test $CXX = cxx; then
      TRILINOS_TEMPLATE_OBJS="cxx_repository/*.o"
    fi
  ;;
esac

AC_MSG_CHECKING([for template-specific object files to be linked])
if test "X$TRILINOS_TEMPLATE_OBJS" = "X"; then
 AC_MSG_RESULT(none)
else
 AC_MSG_RESULT($TRILINOS_TEMPLATE_OBJS)
fi

AC_SUBST(TRILINOS_TEMPLATE_OBJS)

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
##
#
AC_CONFIG_FILES([
		Makefile 
		amd/Makefile
		amd/src/Makefile
		amd/example/Makefile
		src/Makefile
		test/Makefile
		example/Makefile
		example/AmesosDscpack/Makefile
		example/AmesosFactory/Makefile
		example/AmesosSuperludist/Makefile
		test/TestOptions/Makefile
		test/TestWithTriutils/Makefile
		test/scripts/Makefile
		test/scripts/daily/Makefile
		test/scripts/daily/mpi/Makefile
		test/scripts/daily/serial/Makefile
		])

AC_OUTPUT

# Bye World!
echo "----------------------------------------"
echo "Finished Running Amesos Configure Script"
echo "----------------------------------------"
