# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
# 
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2002,2003) Sandia Corporation
# 
#
# Copyright (2001,2002,2003) Sandia Corportation. Under the terms of Contract 
# DE-AC04-94AL85000, there is a non-exclusive license for use of this 
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
# 
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in ths data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
# 
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# Initialization 
# ------------------------------------------------------------------------



# This must be the first line in configure.ac.
# Optional 3rd argument is email address for bugs.

#  #np# - e-mail address below
AC_INIT(new_package, 0.0, YourName@cs.sandia.gov)

# Hello World!
echo "----------------------------------------"
echo "Running New_Package Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.  Any file in that directory will do,
# preferably one that is unlikely to be removed or renamed.

# #np# - file name below
AC_CONFIG_SRCDIR(/src/Newp_Hello.cpp)

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

#
#   I wimped out here.  I still think that ../../config can work 
#   for any package with a frim dependency on Trilinos.
#   But, using ../../config puts some ../../config files into a list
#   of files that it looked like it might delete under certain circumstances,
#   or perhaps put into some tar ball.  Anyway, it scared me off.
#   In order to get rid of the local config directory, you have to 
#   references to the config directory in MAkefile.am and in bootstrap
#   (which should use aclocal -I ../../config)
#
AC_CONFIG_AUX_DIR(config)

# Configure should create src/New_Package_config.h from src/New_Package_config.h.in

AM_CONFIG_HEADER(src/New_Package_config.h:src/New_Package_config.h.in)

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

AM_INIT_AUTOMAKE(new_package, 2.0)

# Specify required version of autoconf.

AC_PREREQ(2.52)

# Specify revision number of configure.ac that was used to generate
# configure.

dnl Commented out this line for now because it causes headaches in
dnl committing.
dnl AC_REVISION($Revision$)

#  #np# - eliminate f77 if your package does not use fortran
# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ gcc c++ cxx)
BASE_CXX=${CXX}
AC_PROG_F77(f77 g77 f90 xlf90)
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI

#   These flags are not fully supported as of Jan 2003
# ------------------------------------------------------------------------
# Checks for --enable-debug, --enable-optimize, --enable-profile, --enable-purify
# Changes flags and compiler names accordingly
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_COMPILER_OPTIONS

# ------------------------------------------------------------------------
# Checks for extra flags
# #np# - You can pull fflags if you aren't using fortran 
# ------------------------------------------------------------------------
# Check if --with-flags present, prepend any specs to FLAGS

TAC_ARG_WITH_FLAGS(ccflags, CCFLAGS)
TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
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
# Special check for case of running on TFLOP
# ------------------------------------------------------------------------
AC_MSG_CHECKING([if host is TFLOP (and therefore needs special TFLOP define statement)])
if test "X${host}" = "Xtflop"; then
  AC_MSG_RESULT(yes)
  AC_DEFINE(TFLOP,,[define if compiling for ASCI Red, aka TFLOP])
else
  AC_MSG_RESULT(no)
fi

#  #np# - May not be required, but strongly recommended.
# ------------------------------------------------------------------------
# Checks for header files - customize Epetra_ConfigDefs.h 
# ------------------------------------------------------------------------

AC_CHECK_HEADERS([cstdlib stdlib.h cstdio stdio.h cassert assert.h],,,)
AC_CHECK_HEADERS([cmath math.h iomanip iomanip.h iostream iostream.h],,,)
AC_CHECK_HEADERS([string string.h],,,)

if test "X$ac_cv_header_cstdlib" != "Xyes" && test "X$ac_cv_header_stdlib_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdlib or stdlib.h.  New_Package requires at least one or the other.]])
fi

if test "X$ac_cv_header_cstdio" != "Xyes" && test "X$ac_cv_header_stdio_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdio or stdio.h.  New_Package requires at least one or the other.]])
fi

if test "X$ac_cv_header_cassert" != "Xyes" && test "X$ac_cv_header_assert_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cassert or assert.h.  New_Package requires at least one or the other.]])
fi

if test "X$ac_cv_header_string" != "Xyes" && test "X$ac_cv_header_string_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find string or string.h.  New_Package requires at least one or the other.]])
fi

if test "X$ac_cv_header_iostream" != "Xyes" && test "X$ac_cv_header_iostream_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find iostream or iostream.h.  New_Package requires at least one or the other.]])
fi

if test "X${host}" != "Xtflop"; then
  if test "X$ac_cv_header_cmath" != "Xyes" && test "X$ac_cv_header_math_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find cmath or math.h.  New_Package requires at least one or the other.]])
  fi
else
  # TFLOP
  if test "X$ac_cv_header_iomanip" != "Xyes" && test "X$ac_cv_header_iomanip_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find iomanip or iomanip.h.  New_Package requires at least one or the other.]])
  fi
fi  

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE

# ------------------------------------------------------------------------
# Checks for compiler characteristics - Epetra_ConfigDefs.h uses these, 
# so it makes sense to leave them in - see PORTABILITY in this directory
# ------------------------------------------------------------------------
AC_CXX_NAMESPACES
AC_CXX_HAVE_STL
AC_CXX_BOOL
AC_CXX_MUTABLE
# ---------------------------------------------------------------------
# These are pretty obscure
# ---------------------------------------------------------------------
AC_CXX_NEW_FOR_SCOPING
#AC_CXX_STD_IOS_BASE_FMTFLAGS
AC_CXX_STD_SPRINTF

# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# #np# newp_swahili is an example of a package feature which requires both
# a cpp #define command and one or more source or header files (and hence
# a change to the makefile).
# If configure is called with "--enable-newp_swahili", it will place 
# a #define HAVE_NEWP_SWAHILI in New_Package_config.h and add source 
# and header files as specified in Makefile.am 
# ------------------------------------------------------------------------
#
TAC_ARG_ENABLE_OPTION(newp_swahili, [enable swahili support], NEWP_SWAHILI, no)
AM_CONDITIONAL(HAVE_NEWP_SWAHILI, [test "X$ac_cv_use_newp_swahili" != "Xno"])

# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------

#  The threads packages do not work yet (Jan 2003)
# TAC_ARG_ENABLE_OPTION(threads, [enable shared memory threads], THREADS, no)
# AM_CONDITIONAL(HAVE_THREADS, [test "X$ac_cv_use_threads" != "Xno"])

#
#  I am not sure what these packages do
TAC_ARG_ENABLE_OPTION(fatal_messages, [enable printing of fatal error messages], 
                   FATAL_MESSAGES, yes)
TAC_ARG_ENABLE_OPTION(warning_messages, [enable printing of warning error messages], 
                   WARNING_MESSAGES, no)
TAC_ARG_ENABLE_OPTION(format_io, [Use package-defined IO formatting for output], 
                   FORMAT_IO, yes)

# These are specific to Epetra
# TAC_ARG_WITH_PACKAGE(zoltan, [Enable Zoltan interface support], ZOLTAN, no)
# TAC_ARG_WITH_PACKAGE(rtop, [Enable RTOp interface support], RTOP, no)
# TAC_ARG_WITH_PACKAGE(esi, [Enable ESI interface support], ESI, yes)
# TAC_ARG_WITH_PACKAGE(superlu, [Enable SuperLU interface support], SUPERLU, no)

# ------------------------------------------------------------------------
# Specify other directories 
# ------------------------------------------------------------------------

# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS

# #np# - Yet another opportunity to remove code if you aren't
# using Fortran
# Define F77_FUNC that will be used to link with Fortran subroutines, e.g., BLAS/LAPACK - trash WORKGXX 
AC_F77_WRAPPERS

# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

# #np# - 
# To add any libraries that you want to check for, replace pow below 
# with the name of an entry point in that library.
# checks for libraries now
dnl Replace `main' with a function in -lm:
AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))

dnl Check for user-specified names for LAPACK and BLAS
TAC_ARG_WITH_BLASLIB
TAC_ARG_WITH_LAPACKLIB


# ------------------------------------------------------------------------
# Checks for compiler characteristics
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for linker characteristics
# ------------------------------------------------------------------------

# Determine libraries needed for linking with Fortran
AC_F77_LIBRARY_LDFLAGS


#  ----------------------------------------------------
#  Needs this because Dec handles templates differently
#  ----------------------------------------------------
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
#  #np# - Most packages should have an example and a full suite 
#  of tests.  New package comes only with src and example.
AC_CONFIG_FILES([
		Makefile 
		src/Makefile
		example/Makefile
		])

AC_OUTPUT()

# Bye World!
echo "----------------------------------------"
echo "Finished Running New_Package Configure Script"
echo "----------------------------------------"
