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

AC_INIT([amesos],[1.1d],[kstanley@cs.utk.edu])

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
TAC_ARG_WITH_LIBS(libs, LIBS)
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

AC_CHECK_HEADERS([valgrind.h valgrind/valgrind.h],,,)

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
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( amesos, tests, [Build Amesos tests], AMESOS_TESTS)
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_amesos_tests" != "Xno")

TAC_ARG_ENABLE_FEATURE(examples, [Build examples for all Trilinos packages (not all packages are sensitive to this option)], EXAMPLES, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( amesos, examples, [Build Amesos examples], AMESOS_EXAMPLES)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_amesos_examples" != "Xno")

TAC_ARG_ENABLE_FEATURE(valgrind, [Use valgrind to perform memory checks and detect leaks (not all packages are sensitive to this option)], VALGRIND, no)
AM_CONDITIONAL(HAVE_VALGRIND, test "X$ac_cv_use_valgrind" != "Xno")

TAC_ARG_ENABLE_FEATURE(libcheck, [Check for some third-party libraries including BLAS and LAPACK.  (Cannot be disabled unless tests and examples are also disabled.)], LIBCHECK, yes)

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

# If tests, examples and libcheck are disabled, we don't have to check
# for these libraries.
if test "X$ac_cv_use_amesos_tests" != "Xno" || test "X$ac_cv_use_amesos_examples" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then

dnl Check for user-specified names for LAPACK and BLAS
TAC_ARG_WITH_BLASLIB
# ------------------------------  ABC
TAC_ARG_WITH_LAPACKLIB
# ------------------------------  ABCD

  ACX_BLAS([], [AC_MSG_ERROR([[Cannot find blas library, specify a path using --with-ldflags=-I<DIR> (ex. --with-ldflags=-I/usr/path/lib) or a specific library using --with-blas=DIR/LIB (for example --with-blas=/usr/path/lib/libcxml.a).  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]])])
# ------------------------------  ABCDE
  ACX_LAPACK([], [AC_MSG_ERROR([[Cannot find lapack library, specify a path using --with-ldflags=-I<DIR> (ex. --with-dlflags=-I/usr/path/lib) or a specific library using --with-lapack=DIR/LIB (ex. --with-lapack=/usr/path/lib/libcxml.a).  See http://software.sandia.gov/trilinos/packages/amesos/ for details.]])])

fi
# end of the list of libraries that don't need to be checked for if
# tests and examples are disabled.

# ------------------------------------------------------------------------
# Checks for addons
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE_SUB( amesos, superludist, [Enable Amesos interface support for SuperLU Dist], 
		AMESOS_SUPERLUDIST, no)

TAC_ARG_ENABLE_FEATURE_SUB( amesos, dscpack, [Enable Amesos interface support for DSCPACK], 
		AMESOS_DSCPACK, no)

TAC_ARG_ENABLE_FEATURE_SUB( amesos, scalapack, [Enable Amesos interface support for SCALAPACK], 
		AMESOS_SCALAPACK, no)

TAC_ARG_ENABLE_FEATURE_SUB( amesos, merikos, [Enable Amesos interface support for MERIKOS], 
		AMESOS_MERIKOS, no)

TAC_ARG_ENABLE_FEATURE_SUB( amesos, umfpack, [Enable Amesos interface support for UMFPACK], 
		AMESOS_UMFPACK, no)

TAC_ARG_ENABLE_FEATURE_SUB( amesos, superlu, [Enable Amesos interface support for SUPERLU], 
		AMESOS_SUPERLU, no)

# always enable klu
TAC_ARG_ENABLE_FEATURE_SUB( amesos, klu, [Enable Amesos interface support for KLU], 
 		AMESOS_KLU, yes)

TAC_ARG_ENABLE_FEATURE_SUB( amesos, mumps, [Enable Amesos interface support for the double-precision version of MUMPS], 
		AMESOS_MUMPS, no)

TAC_ARG_ENABLE_FEATURE_SUB( amesos, smumps, [Enable Amesos interface support for the single-precision version of MUMPS. This is intended for use in conjuction with multilevel preconditioners (in particular ML). Users may decide to use single-precision solves of the coarse problem to save memory and computational time. Note: As Amesos is based on the Epetra_LinearProblem class (defined for double precision only), this interface still requires double-precision matrix and vectors. After the solver phase, the single precision vector is copied into the double-precision solution vector of the given Epetra_LinearProblem. If enabled, disabled double-precision version of MUMPS.], 
		AMESOS_SMUMPS, no)

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
#  Ken Stanley removes this in April 2004 in the hopes that this will allow 
#  amesos to compile and link on Stratus. 
#
#  case "$target" in
#  alpha*-dec-osf*)  
#    if test $CXX = cxx; then
#      CXXFLAGS="-tweak ${CXXFLAGS}"
#    fi
#  ;;
#esac

# removed:
#      CXXFLAGS="  ${CXXFLAGS} -I\${top_srcdir}/amd/src/ -I\${top_srcdir}/src/klu"
#      CFLAGS="  ${CFLAGS} -I\${top_srcdir}/amd/src/"

AC_LANG_PUSH([C])
#  Atlantis fails unless we do the test links with the C compiler (instead of c++)
# AC_LANG([C])

# Fix libraries 

LIBS="${BLAS_LIBS} ${LIBS} ${FLIBS}"
if test "X${HAVE_PKG_MPI}" != "Xyes"; then
  LIBS = "${LIBS} ${MPI_LIBS}"
fi

#  =====  #
#  K L U  #
#  =====  #

if test "X$ac_cv_use_amesos_klu" != "Xno" ; then
  AM_CONDITIONAL(HAVE_AMESOS_KLU, true ) 
  AM_CONDITIONAL(NEED_PARAMETER_LIST, true ) 
else
   AM_CONDITIONAL(HAVE_AMESOS_KLU, false ) 
fi

#  =============  #
#  U M F P A C K  #
#  =============  #

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

#   This causes problems when you use build with both superlu and superludist
#
#    AC_TRY_LINK_FUNC( set_default_options,
#                        AC_DEFINE(OLD_SUPERLUDIST,,[Define if we are using the old Superludist]) ) 

  AC_TRY_LINK_FUNC(umfdi_local_search,,
  AC_MSG_ERROR([[Cannot find the UMFPACK library, specify a location using --with-libs, (ex. --with-libs=/home/username/lib/libumfpack.a, or use --with-libs=-lumfpack --with-ldflags=-L/home/username/lib).]]))

  AC_CHECK_HEADER([umfpack.h],,AC_MSG_ERROR([[Cannot find header files for UMFPACK (umfpack.h), specify a location using --with-incdirs (ex. --with-incdirs=-I/home/username/include).]]))

  AM_CONDITIONAL(HAVE_AMESOS_UMFPACK, true ) 
  AM_CONDITIONAL(NEED_PARAMETER_LIST, true ) 
else
  AM_CONDITIONAL(HAVE_AMESOS_UMFPACK, false ) 
fi

#  =============  #
#  S U P E R L U
#  =============  #

if test "X$ac_cv_use_amesos_superlu" != "Xno"; then

  AC_TRY_LINK_FUNC(SuperLU_timer_,,
  AC_MSG_ERROR([[Cannot find the SuperLU library, specify a location using --with-libs, (ex.  --with-libs=/home/username/lib/libsuperlupack.a, or use --with-libs=-lsuperlupack --with-ldflags=-L/home/username/lib).]]))  

  AC_CHECK_HEADER([Cnames.h],,AC_MSG_ERROR([[Cannot find header files for SuperLU (Cnames.h), specify a location using --with-incdirs (ex.  --with-incdirs=-I/home/username/include).]])) 

  AM_CONDITIONAL(HAVE_AMESOS_SUPERLU, true ) 
  AM_CONDITIONAL(NEED_PARAMETER_LIST, true ) 
else
  AM_CONDITIONAL(HAVE_AMESOS_SUPERLU, false ) 
fi

#  =====================  #
#  S U P E R L U D I S T  #
#  =====================  #

if test "X$ac_cv_use_amesos_superludist" != "Xno"; then

  if test "X${HAVE_PKG_MPI}" != "Xyes"; then
    AC_MSG_ERROR([[Amesos_Superludist requires MPI.]])
  fi

  AC_TRY_LINK_FUNC(pdgssvx,,
  AC_MSG_ERROR([[Cannot find the SuperLU_dist library, specify a location using --with-libs, (ex.  --with-libs=/home/username/lib/libsuperludist.a, or use --with-libs=-lsuperludist --with-ldflags=-L/home/username/lib).]]))  

  AC_CHECK_HEADER([superlu_ddefs.h],,AC_MSG_ERROR([[Cannot find SuperLU_dist header files (superlu_ddefs.h), specify a location using --with-incdirs (ex.  --with-incdirs=/home/username/include).]]))

  AM_CONDITIONAL(HAVE_AMESOS_SUPERLUDIST, true ) 
  AM_CONDITIONAL(NEED_PARAMETER_LIST, true ) 
else
  AM_CONDITIONAL(HAVE_AMESOS_SUPERLUDIST, false ) 
fi

#  =============  #
#  D S C P A C K  #
#  =============  #

if test "X$ac_cv_use_amesos_dscpack" != "Xno"; then

  AC_TRY_LINK_FUNC(DSC_CloseSolver,,
  AC_MSG_ERROR([[Cannot find the DSCPACK library, specify a location using --with-libs, (ex.  --with-libs=/home/username/lib/libdscpack.a, or use --with-libs=-ldscpack --with-ldflags=-L/home/username/lib).]]))  

  AC_CHECK_HEADER([dscmain.h],,AC_MSG_ERROR([[Cannot find DSCPACK header files (dscmain.h), specify a location using --with-incdirs (ex., --with-incdirs=-I/home/username/include).]]))

  AM_CONDITIONAL(HAVE_AMESOS_DSCPACK, true ) 
  AM_CONDITIONAL(NEED_PARAMETER_LIST, true ) 
else
  AM_CONDITIONAL(HAVE_AMESOS_DSCPACK, false ) 
fi

#  =========  #
#  B L A C S  #
#  =========  #

if test "X$ac_cv_use_amesos_scalapack" != "Xno" || 
   test "X$ac_cv_use_amesos_mumps" != "Xno"; then

  AC_F77_FUNC(dgebs2d)

  AC_TRY_LINK_FUNC($dgebs2d,,
  AC_MSG_ERROR([[Cannot find the BLACS library (used by MUMPS and/or ScaLAPACK), specify a location using for example --with-libs=/home/username/lib/libblacs.a, or --with-libs=-lblacs --with-ldflags=-L/home/username/lib).]]))  

fi
		
#  ===========  #
#  C B L A C S  #
#  ===========  #

if test "X$ac_cv_use_amesos_scalapack" != "Xno" || 
   test "X$ac_cv_use_amesos_mumps" != "Xno"; then

  AC_TRY_LINK_FUNC(Cdgebs2d,,
  AC_MSG_ERROR([[Cannot find the CBLACS library (used by MUMPS and/or ScaLAPACK), specify a location using for example --with-libs=/home/username/lib/libblacs.a, or --with-libs=-lblacs --with-ldflags=-L/home/username/lib).]]))  

fi
		
#  =================  #
#  S C A L A P A C K  #
#  =================  #

if test "X$ac_cv_use_amesos_scalapack" != "Xno"; then

  AC_F77_FUNC(pdgetrf)

  AC_TRY_LINK_FUNC($pdgetrf,,
  AC_MSG_ERROR([[Cannot find the SCALAPACK library, specify a location using --with-libs, (ex.  --with-libs=/home/username/lib/libscalack.a, or use --with-libs=-lscalapack --with-ldflags=-L/hhome/username/lib).]]))  

  AM_CONDITIONAL(HAVE_AMESOS_SCALAPACK, true ) 
  AM_CONDITIONAL(NEED_PARAMETER_LIST, true ) 
else
  AM_CONDITIONAL(HAVE_AMESOS_SCALAPACK, false ) 
fi

#  =========  #
#  M U M P S  #
#  =========  #

# PORD (required by both MUMPS and SMUMPS)

if test "X$ac_cv_use_amesos_mumps" != "Xno" || 
  test "X$ac_cv_use_amesos_smumps" != "Xno" ; then

  AC_TRY_LINK_FUNC(newMultisector,,
  AC_MSG_ERROR([[Cannot find the PORD library (the default ordering, required by MUMPS), specify a location using --with-libs, (ex.  --with-libs=/home/username/lib/libpord.a, or use --with-libs=-lpord --with-ldflags=-L/home/username/lib).]]))  

fi

# DMUMPS

if test "X$ac_cv_use_amesos_mumps" != "Xno"; then

  AC_TRY_LINK_FUNC(dmumps_c,, AC_MSG_ERROR([[Cannot find the double precision MUMPS library, specify a location using --with-libs, (ex.  --with-libs=/home/username/lib/libdmumps.a, or use --with-libs=-ldmumps --with-ldflags=-L/hhome/username/lib).]]))  

AC_CHECK_HEADER([dmumps_c.h],,AC_MSG_ERROR([[Cannot find header files for MUMPS (dmumps_c.h), specify a location using --with-incdirs (ex.  --with-incdirs=-I/home/username/include).]]))

  AM_CONDITIONAL(HAVE_AMESOS_MUMPS, true ) 
  AM_CONDITIONAL(NEED_PARAMETER_LIST, true ) 
else
  AM_CONDITIONAL(HAVE_AMESOS_MUMPS, false ) 
fi

# SMUMPS

if test "X$ac_cv_use_amesos_smumps" != "Xno"; then

  AC_TRY_LINK_FUNC(smumps_c,,
  AC_MSG_ERROR([[Cannot find the single precision MUMPS library, specify a location using --with-libs, (ex.  --with-libs=/home/username/lib/libsmumps.a, or use --with-libs=-lsmumps --with-ldflags=-L/home/username/lib).]]))  

  AC_CHECK_HEADER([smumps_c.h],,AC_MSG_ERROR([[Cannot find header files for MUMPS (smumps_c.h), specify a location using --with-incdirs (ex.  --with-incdirs=-I/home/username/include).]]))

  AM_CONDITIONAL(HAVE_AMESOS_SMUMPS, true ) 
  AM_CONDITIONAL(NEED_PARAMETER_LIST, true ) 
else
  AM_CONDITIONAL(HAVE_AMESOS_SMUMPS, false ) 
fi

#  =============  #
#  M E R I K O S  #
#  =============  #

if test "X$ac_cv_use_amesos_merikos" != "Xno"; then
  AM_CONDITIONAL(HAVE_AMESOS_MERIKOS, true ) 
  AM_CONDITIONAL(NEED_PARAMETER_LIST, true ) 
else
  AM_CONDITIONAL(HAVE_AMESOS_MERIKOS, false ) 
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


AC_LANG_POP

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
