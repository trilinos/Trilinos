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


#
#  Changes made April 2004:
#    Each Amesos class will:
#       1)  First try to compile without specifying the library.  This allows
#           the user to specify the needed third party library at compile time.
#       2)  Next try to compile with -lsuperludist.  This allows the user to specify 
#           the directory while using the default name for the library
#       3)  Then try the existing Triliinos3PL scheme. 
#

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

TAC_ARG_WITH( trilinos3pldir, [Set Trilinos3PL directory], trilinos3pldir_not_specified )

## echo "  tac_with_trilinos3pldir =  ${tac_with_trilinos3pldir} " 
TAC_CHECK_FILE(${tac_with_trilinos3pldir})

#
#  --enable-amesos-slud, --with-amesos-sluddir, --with-amesos-sludlib, 
#    --with-amesos-sludincdir
#
TAC_ARG_ENABLE_FEATURE_SUB( amesos, superludist, [Enable Amesos interface support for SuperLU Dist], 
		AMESOS_SUPERLUDIST, no)

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
    AC_MSG_ERROR([[Amesos_Superludist requires MPI.  See http://software.sandia.gov/trilinos/packages/amesos/ for details]])
  fi

  TAC_ARG_WITH_SUB( amesos, superludistdir, [Set SuperLU dist directory],  ${tac_with_trilinos3pldir}/SuperLU_DIST_2.0 )
  TAC_CHECK_FILE(${tac_with_amesos_superludistdir})

#
#       1)  First try to compile without specifying the library.  This allows
#           the user to specify the needed third party library at compile time.
#
  SAVE_LIBS="${LIBS}"
  LIBS="${BLAS_LIBS} ${MPI_LIBS} ${FLIBS} ${LIBS}"
  AC_TRY_LINK_FUNC(pdgssvx,superludist_lib_ok="yes",superludist_lib_ok="no")
     


  if test "X$superludist_lib_ok" != "Xyes"; then
#
#       2)  Next try to compile with -lsuperludist.  This allows the user to specify 
#           the directory while using the default name for the library
     LIBS=" -lsuperludist ${LIBS}"

     AC_TRY_LINK_FUNC(pdgssvx,superludist_with_lib="yes",superludist_with_lib="no")

     if test "X$superludist_with_lib" == "Xyes"; then
       amesos_libs="${amesos_libs} -lsuperludist"
     else
#       3)  Then try the existing Triliinos3PL scheme. 
#
       TAC_ARG_WITH_SUB( amesos, superludistlib, [Set SuperLU dist library],  
                         ${tac_with_amesos_superludistdir}/libsuperludist.a )
       TAC_CHECK_FILE(${tac_with_amesos_superludistlib})
       LIBS="${tac_with_amesos_superludistlib} ${BLAS_LIBS} ${MPI_LIBS} ${FLIBS} ${SAVE_LIBS}"
       amesos_libs="${amesos_libs} ${tac_with_amesos_superludistlib}"
     fi
  fi
#
#   At this point, LIBS and amesos_libs both have our best shot, if they don't work now
#   we give up
#
   AC_TRY_LINK_FUNC( set_default_options,
                       AC_DEFINE(OLD_SUPERLUDIST,,[Define if we are using the old Superludist]) ) 

   AC_TRY_LINK_FUNC(pdgssvx,,
      AC_MSG_ERROR([[Cannot find SuperLU Dist library, specify a location using --with-amesos-superludistlib, (ex. --with-amesos-superludistlib=/home/username/Trilinos3PL/SuperLU_DIST/libsuperludist.a)]]))  
    LIBS="${SAVE_LIBS}"

#
#  First try to find superlu_ddefs.h without using the SuperLU dist incude directory
#
  AC_CHECK_HEADERS([superlu_ddefs.h],,,)

  if test "X$ac_cv_header_superlu_ddefs_h" != "Xyes"; then
#
#   if we can't find supelu_ddefs.h in the include search path, look for it in superludistincir
#   which defaults to Trilinos3PL/SuperLU_DIST_2.0/SRC
#
    TAC_ARG_WITH_SUB( amesos, superludistincdir, [Set SuperLU dist include directory],  ${tac_with_amesos_superludistdir}/SRC )

    TAC_CHECK_FILE(${tac_with_amesos_superludistincdir}/superlu_ddefs.h,,
      AC_MSG_ERROR([[Cannot find superlu_ddefs.h in ${tac_with_amesos_superludistincdir} , specify a location using --with-amesos-superludistincdir, --with-amesos-superludistdir or --with-trilinos3pldir.  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]]) )

    CXXFLAGS="-I${tac_with_amesos_superludistincdir} ${CXXFLAGS}"
  fi

  AM_CONDITIONAL(HAVE_AMESOS_SUPERLUDIST, true ) 
  AM_CONDITIONAL(NEED_PARAMETER_LIST, true ) 
else
  AM_CONDITIONAL(HAVE_AMESOS_SUPERLUDIST, false ) 
fi



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

#
#       1)  First try to compile without specifying the library.  This allows
#           the user to specify the needed third party library at compile time.
#
  SAVE_LIBS="${LIBS}"
  LIBS="${BLAS_LIBS} ${MPI_LIBS} ${FLIBS} ${LIBS}"
  AC_TRY_LINK_FUNC(DSC_CloseSolver,dscpack_lib_ok="yes",dscpack_lib_ok="no")
     AC_MSG_CHECKING( 465 ) 
     AC_MSG_RESULT( 465 ) 
     
  if test "X$dscpack_lib_ok" != "Xyes"; then
#
#       2)  Next try to compile with -ldscdbl.  This allows the user to specify 
#           the directory while using the default name for the library
     LIBS=" -ldscdbl ${LIBS}"

     AC_MSG_CHECKING( 472 ) 
     AC_MSG_RESULT( 472 ) 
     AC_TRY_LINK_FUNC(DSC_CloseSolver,dscpack_lib_ok="yes",dscpack_lib_ok="no")

     if test "X$dscpack_lib_ok" == "Xyes"; then
     AC_MSG_CHECKING( 480 ) 
     AC_MSG_RESULT( 480 ) 
       amesos_libs="${amesos_libs} -ldscdbl"
     else
#       3)  Then try the existing Triliinos3PL scheme with dsclibdbl.a as the default library
#
       TAC_ARG_WITH_SUB( amesos, dscpacklib, [Set DSCPACK library],  
                         ${tac_with_amesos_dscpackdir}/DSC_LIB/dsclibdbl.a )
       TAC_CHECK_FILE(${tac_with_amesos_dscpacklib})
       LIBS="${tac_with_amesos_dscpacklib} ${BLAS_LIBS} ${MPI_LIBS} ${FLIBS} ${SAVE_LIBS}"
       AC_TRY_LINK_FUNC(DSC_CloseSolver,dscpack_lib_ok="yes",dscpack_lib_ok="no")
        if test "X$dscpack_lib_ok" == "Xyes"; then
     AC_MSG_CHECKING( 488 ) 
     AC_MSG_RESULT( 488 ) 
          amesos_libs="${amesos_libs} ${tac_with_amesos_dscpacklib}"
       else
#         4)  Then try  the existing Triliinos3PL scheme with libdscdbl.a as the default library
#
     AC_MSG_CHECKING( 493 ) 
     AC_MSG_RESULT( 493 ) 
          LIBS="${tac_with_amesos_dscpackdir}/DSC_LIB/libdscdbl.a ${BLAS_LIBS} ${MPI_LIBS} ${FLIBS} ${SAVE_LIBS}"
          amesos_libs="${amesos_libs} ${tac_with_amesos_dscpackdir}/DSC_LIB/libdscdbl.a"
       fi
     fi
  fi
#
#   At this point, LIBS and amesos_libs both have our best shot, if they don't work now
#   we give up
#

     AC_MSG_RESULT( " :::509 LIBS = ${LIBS} " ) 

     AC_TRY_LINK_FUNC(DSC_CloseSolver,,
      AC_MSG_ERROR([[Cannot find the DSCPACK library, specify a location using --with-amesos-dscpacklib, (ex. --with-amesos-dscpacklib=/home/username/Trilinos3PL/DSCPACK1.0/DSC_LIB/dsclibdbl.a).  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]]))  
    LIBS="${SAVE_LIBS}"

#
#  First try to find dscmain.h without using the DSCPACK incude directory
#
  AC_CHECK_HEADERS([dscmain.h],,,)

  if test "X$ac_cv_header_dscmain_h" != "Xyes"; then
#
#   if we can't find dscmain.h in the include search path, look for it in dscmainincdir
#   which defaults to Trilinos3PL/DSCPACK1.0/SRC
#
    TAC_ARG_WITH_SUB( amesos, dscpackincdir, [Set DSCPACK include directory],  ${tac_with_amesos_dscpackdir}/DSC_LIB )

    TAC_CHECK_FILE(${tac_with_amesos_dscpackincdir}/dscmain.h,,
      AC_MSG_ERROR([[Cannot find dscmain.h in ${tac_with_amesos_dscpackincdir} , specify a location using --with-amesos-dscpackincdir, --with-amesos-dscpackdir or --with-trilinos3pldir.  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]]) )

    CXXFLAGS="-I${tac_with_amesos_dscpackincdir} ${CXXFLAGS}"
  fi

  AM_CONDITIONAL(HAVE_AMESOS_DSCPACK, true ) 
  AM_CONDITIONAL(NEED_PARAMETER_LIST, true ) 
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

#
#       1)  First try to compile without specifying the library.  This allows
#           the user to specify the needed third party library at compile time.
#
  SAVE_LIBS="${LIBS}"
  LIBS="${BLAS_LIBS} ${MPI_LIBS} ${FLIBS} ${LIBS}"
  AC_TRY_LINK_FUNC(umfdi_local_search,umfpack_lib_ok="yes",umfpack_lib_ok="no")
     
  AC_MSG_RESULT( checkpoint umfpack  A ) 
  AC_MSG_RESULT( checkpoint umfpack  A ) 
  if test "X$umfpack_lib_ok" != "Xyes"; then
#
#       2)  Next try to compile with -lumfpack.  This allows the user to specify 
#           the directory while using the default name for the library
     LIBS=" -lumfpack ${LIBS}"

     AC_MSG_RESULT( checkpoint umfpack  B ) 
     AC_MSG_RESULT( checkpoint umfpack  B ) 
     AC_TRY_LINK_FUNC(umfdi_local_search,umfpack_with_lib="yes",umfpack_with_lib="no")

     if test "X$umfpack_with_lib" == "Xyes"; then
       amesos_libs="${amesos_libs} -lumfpack"
       AC_MSG_RESULT( checkpoint umfpack  C ) 
       AC_MSG_RESULT( checkpoint umfpack  C ) 
     else
#       3)  Then try the existing Triliinos3PL scheme. 
#
       AC_MSG_RESULT( checkpoint umfpack  D ) 
       AC_MSG_RESULT( checkpoint umfpack  D ) 
       TAC_ARG_WITH_SUB( amesos, umfpacklib, [Set UMFPACK library],  
                         ${tac_with_amesos_umfpackdir}/UMFPACK/Lib/libumfpack.a )
       TAC_CHECK_FILE(${tac_with_amesos_umfpacklib})
       LIBS="${tac_with_amesos_umfpacklib} ${BLAS_LIBS} ${MPI_LIBS} ${FLIBS} ${SAVE_LIBS}"
       amesos_libs="${amesos_libs} ${tac_with_amesos_umfpacklib}"
     fi
  fi

  AC_TRY_LINK_FUNC(umfdi_local_search,,
  AC_MSG_ERROR([[Cannot find the UMFPACK library, specify a location using --with-amesos-umfpacklib, (ex. --with-amesos-umfpacklib=/home/username/Trilinos3PL/UMFPACKv4.1/libumfpack.a).  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]]))  



#
#       1)  First try to compile without specifying the library.  This allows
#           the user to specify the needed third party library at compile time.
#
  SAVE_LIBS="${LIBS}"
  AC_TRY_LINK_FUNC(amd_defaults,amd_lib_ok="yes",amd_lib_ok="no")
     
  AC_MSG_RESULT( checkpoint amd  A ) 
  AC_MSG_RESULT( checkpoint amd  A ) 
  if test "X$amd_lib_ok" != "Xyes"; then
#
#       2)  Next try to compile with -lamd.  This allows the user to specify 
#           the directory while using the default name for the library
     LIBS=" -lamd ${LIBS}"

     AC_MSG_RESULT( checkpoint amd  B ) 
     AC_MSG_RESULT( checkpoint amd  B ) 
     AC_TRY_LINK_FUNC(amd_defaults,amd_with_lib="yes",amd_with_lib="no")

     if test "X$amd_with_lib" == "Xyes"; then
       AC_MSG_RESULT( checkpoint amd  C ) 
       AC_MSG_RESULT( checkpoint amd  C ) 
       amesos_libs="${amesos_libs} -lamd"
     else
#       3)  Then try the existing Triliinos3PL scheme. 
#
       AC_MSG_RESULT( checkpoint amd D ) 
       AC_MSG_RESULT( checkpoint amd D ) 
       TAC_ARG_WITH_SUB( amesos, umfpackamdlib, [Set AMD library],  
                         ${tac_with_amesos_umfpackdir}/AMD/Lib/libamd.a )
       TAC_CHECK_FILE(${tac_with_amesos_umfpackamdlib})
       LIBS="${tac_with_amesos_umfpackamdlib} ${BLAS_LIBS} ${MPI_LIBS} ${FLIBS} ${SAVE_LIBS}"
       amesos_libs="${amesos_libs} ${tac_with_amesos_umfpackamdlib}"
     fi
  fi

  AC_TRY_LINK_FUNC(amd_defaults,,
  AC_MSG_ERROR([[Cannot find UMFPACK's AMD library, specify a location using --with-amesos-umfpackamdlib, (ex. --with-amesos-umfpackamdlib=/home/username/Trilinos3PL/AMD/Lib/libamd.a).  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]]))  

  LIBS="${SAVE_LIBS}"

#
#  First try to find umfpack.h without using the UMFPACK incude directory
#
  AC_CHECK_HEADERS([umfpack.h],,,)

  if test "X$ac_cv_header_umfpack_h" != "Xyes"; then
#
#   if we can't find umfpack.h in the include search path, look for it in umfpackincir
#   which defaults to Trilinos3PL/UMFPACKv4.1/UMFPACK/Include
#
    TAC_ARG_WITH_SUB( amesos, umfpackincdir, [Set UMFPACK include directory],  
                      ${tac_with_amesos_umfpackdir}/UMFPACK/Include )

    TAC_CHECK_FILE(${tac_with_amesos_umfpackincdir}/umfpack.h,,
      AC_MSG_ERROR([[Cannot find umfpack.h in ${tac_with_amesos_umfpackincdir} , specify a location using --with-amesos-umfpackincdir, --with-amesos-umfpackdir or --with-trilinos3pldir.  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]]) )

    CXXFLAGS="-I${tac_with_amesos_umfpackincdir} ${CXXFLAGS}"
  fi


#
#  First try to find amd.h without using the UMFPACK AMD incude directory
#
  AC_CHECK_HEADERS([amd.h],,,)

  if test "X$ac_cv_header_amd_h" != "Xyes"; then
#
#   if we can't find amd.h in the include search path, look for it in umfpackamdincdir
#   which defaults to Trilinos3PL/UMFPACKv4.1/AMD/Include
#
    TAC_ARG_WITH_SUB( amesos, umfpackamdincdir, [Set UMFPACK include directory],  ${tac_with_amesos_umfpackdir}/AMD/Include )

    TAC_CHECK_FILE(${tac_with_amesos_umfpackamdincdir}/amd.h,,
      AC_MSG_ERROR([[Cannot find amd.h in ${tac_with_amesos_umfpackamdincdir} , specify a location using --with-amesos-umfpackamdincdir, --with-amesos-umfpackdir or --with-trilinos3pldir.  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]]) )

    CXXFLAGS="-I${tac_with_amesos_umfpackamdincdir} ${CXXFLAGS}"
  fi

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



#   April 23, 2004 - None of the blacs, cblacs, scalapack or mumps stuff is tested at this point.



#
#  Look for scalapack, blacs and cblacs libraries 
#    if --enable-amesos-scalapack or --enable-amesos-mumps
#
if test "X$ac_cv_use_amesos_scalapack" != "Xno" || 
   test "X$ac_cv_use_amesos_mumps" != "Xno"; then

#
#   BLACS library
#
  TAC_ARG_WITH_SUB( amesos, blacsdir, [Set BLACS directory],  ${tac_with_trilinos3pldir}/BLACS/LIB )
  TAC_CHECK_FILE(${tac_with_amesos_blacsdir})

  SAVE_LIBS="${LIBS}"
  LIBS="${BLAS_LIBS} ${MPI_LIBS} ${FLIBS} ${LIBS}"
  BASELIBS="${LIBS}"
  BLACSLIB=""
  AC_F77_FUNC(dgebs2d)

  AC_TRY_LINK_FUNC($dgebs2d,blacs_lib_ok="yes",blacs_lib_ok="no")
     
  AC_MSG_RESULT( checkpoint blacs  A ) 
  AC_MSG_RESULT( checkpoint blacs  A ) 
  if test "X$blacs_lib_ok" != "Xyes"; then
     LIBS=" -lblacs ${BASELIBS}"

     AC_MSG_RESULT( checkpoint blacs  B ) 
     AC_MSG_RESULT( checkpoint blacs  B ) 
     AC_TRY_LINK_FUNC($dgebs2d,blacs_with_lib="yes",blacs_with_lib="no")

     if test "X$blacs_with_lib" == "Xyes"; then
       BLACSLIB=" -lblacs "
       AC_MSG_RESULT( checkpoint blacs  C ) 
       AC_MSG_RESULT( checkpoint blacs  C ) 
     else
       AC_MSG_RESULT( checkpoint blacs  D ) 
       AC_MSG_RESULT( checkpoint blacs  D ) 
       TAC_ARG_WITH_SUB( amesos, blacslib, [Set BLACS library],  
                         ${tac_with_amesos_blacsdir}/libblacs.a )
       LIBS=" ${tac_with_amesos_blacslib} ${BASELIBS}"
       AC_TRY_LINK_FUNC($dgebs2d, BLACSLIB=${tac_with_amesos_blacslib}, ) 
       AC_MSG_RESULT( "checkpoint blacs  E BLACSLIB= ${tac_with_amesos_blacslib} " ) 
     fi
  fi

       AC_MSG_RESULT( "checkpoint yo amesos_libs = ${amesos_libs}" ) 
       AC_MSG_RESULT( "checkpoint yo BLACSLIB = ${BLACSLIB}" ) 

  LIBS="${BLACSLIB} ${BASELIBS}"
  AC_TRY_LINK_FUNC($dgebs2d,,
  AC_MSG_ERROR([[Cannot find the BLACS library, specify a location using --with-amesos-blacslib, (ex. --with-amesos-blacslib=/home/username/Trilinos3PL/BLACS/LIB/libblacs.a).  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]]))  


#
#      CBLACS library
#  Look for the cblacs library - some implementations of either ScaLAPACK or MUMPS may need 
#  the cblacs library.  Other implementations do not.  If we find the cblacs, we include it 
#  in the search path, if we don't, we hope for the best.
#
  AC_TRY_LINK_FUNC(Cblacs_pinfo,cblacs_lib_ok="yes",cblacs_lib_ok="no")
  AC_MSG_RESULT( checkpoint cblacs  A ) 
  AC_MSG_RESULT( checkpoint cblacs  A ) 
  if test "X$cblacs_lib_ok" != "Xyes"; then
     LIBS=" -lcblacs ${BLACSLIB} ${BASELIBS}"

     AC_MSG_RESULT( checkpoint cblacs  B ) 
     AC_MSG_RESULT( checkpoint cblacs  B ) 
     AC_TRY_LINK_FUNC(Cblacs_pinfo,cblacs_with_lib="yes",cblacs_with_lib="no")

     if test "X$cblacs_with_lib" == "Xyes"; then
       CBLACSLIB=" -lcblacs ${BLACSLIB} " 
       AC_MSG_RESULT( "checkpoint yy CBLACSLIB = ${CBLACSLIB}" ) 
       AC_MSG_RESULT( checkpoint cblacs  C ) 
       AC_MSG_RESULT( checkpoint cblacs  C ) 
     else
       AC_MSG_RESULT( checkpoint cblacs  D ) 
       AC_MSG_RESULT( checkpoint cblacs  D ) 
       TAC_ARG_WITH_SUB( amesos, cblacslib, [Set CBLACS library],  
                         ${tac_with_amesos_blacsdir}/libcblacs.a )
       LIBS=" ${tac_with_amesos_cblacslib} ${BLACSLIB} ${BASELIBS}"
       AC_TRY_LINK_FUNC(Cblacs_pinfo, CBLACSLIB=" ${tac_with_amesos_cblacslib}  ${BLACSLIB} ", ) 
       AC_MSG_RESULT( "checkpoint cblacs  E CBLACSLIB=${tac_with_amesos_cblacslib}" ) 
     fi
  fi
fi
#
#  No check here because not all implementations require the cblacs library
#

       AC_MSG_RESULT( "checkpoint yoyo CBLACSLIB = ${CBLACSLIB}" ) 
#
#    SCALAPACK library
#
TAC_ARG_WITH_SUB( amesos, scalapackdir, [Set SCALAPACK directory],  ${tac_with_trilinos3pldir}/SCALAPACK )
if test "X$ac_cv_use_amesos_scalapack" != "Xno"; then

  AC_F77_FUNC(pdgetrf)
  LIBS=" ${BLACSLIB} ${CBLACSLIB} ${BASELIBS}"
  SCALAPACKLIB=""
  AC_TRY_LINK_FUNC($pdgetrf,scalapack_lib_ok="yes",scalapack_lib_ok="no")
  AC_MSG_RESULT( checkpoint scalapack  A ) 
  AC_MSG_RESULT( checkpoint scalapack  A ) 
  if test "X$scalapack_lib_ok" != "Xyes"; then
     LIBS=" -lscalapack ${BLACSLIB} ${CBLACSLIB} ${BASELIBS}"

     AC_MSG_RESULT( checkpoint scalapack  B ) 
     AC_MSG_RESULT( checkpoint scalapack  B ) 
     AC_TRY_LINK_FUNC($pdgetrf,scalapack_with_lib="yes",scalapack_with_lib="no")
     if test "X$scalapack_with_lib" == "Xyes"; then
       SCALAPACKLIB=" -lscalapack "
       AC_MSG_RESULT( checkpoint scalapack  C ) 
       AC_MSG_RESULT( checkpoint scalapack  C ) 
     else
       AC_MSG_RESULT( checkpoint scalapack  D ) 
       AC_MSG_RESULT( checkpoint scalapack  D ) 
       TAC_ARG_WITH_SUB( amesos, scalapacklib, [Set SCALAPACK library],  
                         ${tac_with_amesos_scalapackdir}/libscalapack.a )
       LIBS=" ${tac_with_amesos_scalapacklib} ${CBLACSLIB} ${BLACSLIB} ${BASELIBS}"
       AC_TRY_LINK_FUNC($pdgetrf, SCALAPACKLIB=" ${tac_with_amesos_scalapacklib} ", ) 
       AC_MSG_RESULT( "checkpoint scalapack  E SCALACAPKLIB=${tac_with_amesos_scalapcklib}" ) 
     fi
  fi

  LIBS="${SCALAPACKLIB} ${BLACSLIB} ${CBLACSLIB} ${BASELIBS}"
  AC_MSG_RESULT( "checkpoint yoyoyo CBLACSLIB = ${CBLACSLIB}" ) 
  amesos_libs=" ${SCALAPACKLIB} ${BLACSLIB} ${CBLACSLIB} ${amesos_libs} "
  AC_MSG_RESULT( "LIBS = ${LIBS}" ) 

  AC_TRY_LINK_FUNC($pdgetrf,,
  AC_MSG_ERROR([[Cannot find scalapack library, specify a location using --with-amesos-scalapacklib, (ex. --with-amesos-scalapacklib=/home/username/Trilinos3PL/SuperLU_DIST_2.0/libsuperludist.a).  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]]))  

  LIBS="${SAVE_LIBS}"
  AM_CONDITIONAL(HAVE_AMESOS_SCALAPACK, true ) 
  AM_CONDITIONAL(NEED_PARAMETER_LIST, true ) 

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

  LIBS="${SCALAPACKLIB} ${BLACSLIB} ${CBLACSLIB} ${BASELIBS}"

#
#  We don't check do as thorough a check for Mumps at this point.  In fact, we don't do much of a check 
#  at all.
#
  TAC_ARG_WITH_SUB( amesos, mumpsdir, [Set MUMPS directory],  ${tac_with_trilinos3pldir}/MUMPS_4.2_beta/ )
  TAC_CHECK_FILE(${tac_with_amesos_mumpsdir})

  TAC_ARG_WITH_SUB( amesos, mumpslib, [Set MUMPS library],  ${tac_with_amesos_mumpsdir}/lib/libdmumps.a )
  TAC_CHECK_FILE(${tac_with_amesos_mumpslib})

  TAC_ARG_WITH_SUB( amesos, mumpspordlib, [Set MUMPS PORD library],  ${tac_with_amesos_mumpsdir}/lib/libpord.a )
  TAC_CHECK_FILE(${tac_with_amesos_mumpspordlib})

  TAC_ARG_WITH_SUB( amesos, mumpsincdir, [Set MUMPS include directory],  ${tac_with_amesos_mumpsdir}/include )
  TAC_CHECK_FILE(${tac_with_amesos_mumpsincdir})

#
# This is a bad way to check for dmumps_c.h, it won't pick dmumps_c.h up if
# it is in the include path, i.e. /usr/include - see tac_arg_check_mpi for a better way
#

#
#  First try to find amd.h without using the UMFPACK AMD incude directory
#
  AC_CHECK_HEADERS([dmumps_c.h],,,)

  if test "X$ac_cv_header_dmumps_c_h" != "Xyes"; then
#
#   if we can't find dmumps_c.h in the include search path, look for it in mumpsincdir
#   which defaults to Trilinos3PL/MUMPS_4.2_beta/include
#
  TAC_CHECK_FILE(${tac_with_amesos_mumpsincdir}/dmumps_c.h,,
    AC_MSG_ERROR([[Cannot find dmumps_c.h in ${tac_with_amesos_mumps} , specify a location using --with-amesos-mumpsincdir, --with-amesos-mumpsdir or --with-trilinos3pldir.  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]] ))

    CXXFLAGS="-I${tac_with_amesos_mumpsincdir} ${CXXFLAGS}"
  fi

#
#    The following attempt to confirm that libmumps.a and lidpord.a are the right libraries
#    fails - presumably because we need to link other libraries in (either fortran 
#    scalapack or blacs libraries.)
#

#    SAVE_LIBS="${LIBS}"
#      LIBS="${tac_with_amesos_mumpslib} ${MPI_LIBS} ${LIBS}"
#      AC_TRY_LINK_FUNC(dmumps_164_,,
#      AC_MSG_ERROR([[Cannot find the MUMPS library, specify a location using --with-amesos-mumpslib, (ex. --with-amesos-mumpslib=/home/username/Trilinos3PL/MUMPS_4.2_beta/lib/libdmumps.a).  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]]))  
#    LIBS="${tac_with_amesos_mumpspordlib} ${MPI_LIBS} ${LIBS}"
#      AC_TRY_LINK_FUNC(qsortUpInts,,
#      AC_MSG_ERROR([[Cannot find the PORD library needed by MUMPS, specify a location using --with-amesos-mumpspordlib, (ex. --with-amesos-mumpslib=/home/username/Trilinos3PL/MUMPS_4.2_beta/lib/libpord.a).  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]]))  
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
      AC_MSG_ERROR([[Cannot find SuperLU serial library, specify a location using --with-amesos-superlulib, (ex. --with-amesos-superlulib=/home/username/Trilinos3PL/SuperLU/libsuperlu.a).  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]]))  
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
if test "X$TRILINOS_TEMPLATE_OBJS" = "X"; then
  case $target in
    alpha*)
      if test $CXX = cxx; then
        TRILINOS_TEMPLATE_OBJS="cxx_repository/*.o"
      fi
    ;;
  esac
fi

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
