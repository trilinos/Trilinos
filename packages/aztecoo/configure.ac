# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
# 
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2002) Sandia Corporation
# 
#
# Copyright (2001) Sandia Corportation. Under the terms of Contract 
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

AC_INIT(aztecoo, 2.0, mheroux@cs.sandia.gov)

# Hello World!
echo "----------------------------------------"
echo "Running AztecOO Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.

AC_CONFIG_SRCDIR(/src/az_ilut.c)

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

AC_CONFIG_AUX_DIR(config)

# Configure should create src/AztecOO_config.h from src/AztecOO_config.h.in

AM_CONFIG_HEADER(src/AztecOO_config.h:src/AztecOO_config.h.in)

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

AM_INIT_AUTOMAKE(aztecoo, 2.0)

# Specify required version of autoconf.

AC_PREREQ(2.52)

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
# Checks for programs.
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
# Checks for header files
# ------------------------------------------------------------------------

AC_HEADER_STDC
AC_CHECK_HEADERS([cstdlib stdlib.h cstdio stdio.h cassert assert.h],,,)
AC_CHECK_HEADERS([cmath math.h string string.h cfloat float.h],,,)
AC_CHECK_HEADERS([malloc.h ctime sys/time.h])

if test "X$ac_cv_header_cstdlib" != "Xyes" && test "X$ac_cv_header_stdlib_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdlib or stdlib.h.  AztecOO requires at least one or the other.]])
fi

if test "X$ac_cv_header_cstdio" != "Xyes" && test "X$ac_cv_header_stdio_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdio or stdio.h.  AztecOO requires at least one or the other.]])
fi

if test "X$ac_cv_header_cassert" != "Xyes" && test "X$ac_cv_header_assert_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cassert or assert.h.  AztecOO requires at least one or the other.]])
fi

if test "X$ac_cv_header_string" != "Xyes" && test "X$ac_cv_header_string_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find string or string.h.  AztecOO requires at least one or the other.]])
fi

  if test "X$ac_cv_header_cmath" != "Xyes" && test "X$ac_cv_header_math_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find cmath or math.h.  AztecOO requires at least one or the other.]])
  fi

  if test "X$ac_cv_header_cfloat" != "Xyes" && test "X$ac_cv_header_float_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find cfloat or float.h.  AztecOO requires at least one or the other.]])
  fi

  if test "X$ac_cv_header_ctime" != "Xyes" && test "X$ac_cv_header_sys_time_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find ctime or sys/time.h.  AztecOO requires at least one or the other.]])
  fi

  if test "X$ac_cv_header_malloc_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find malloc.h.  AztecOO requires this header file.]])
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
#AC_CXX_STD_IOS_BASE_FMTFLAGS
AC_CXX_STD_SPRINTF

# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE(oldifpack, [Enable old ifpack interface support], OLDIFPACK, no)
AM_CONDITIONAL(BUILD_OLDIFPACK, test "X$ac_cv_use_oldifpack" != "Xno")

# The following section is required to decide if the AZ_lu preconditioner
# should be built.  If not, y12m is not necessary.
AC_ARG_ENABLE(aztecoo-azlu,
 [AC_HELP_STRING([--enable-aztecoo-azlu],[Enable az-lu preconditioner.  Default is no.  Requires y12m.])],
 ac_cv_use_aztecoo_azlu=$enableval, ac_cv_use_aztecoo_azlu=no)

AC_MSG_CHECKING(whether to use aztecoo-azlu)

if test "X$ac_cv_use_aztecoo_azlu" != "Xno"; then
  AC_MSG_RESULT(yes)
  AC_DEFINE([HAVE_AZLU],,[Define if you want to build AZ_lu])
else
  AC_MSG_RESULT(no)
fi

dnl THIS NEEDS TO BE FILLED IN!!!

# ------------------------------------------------------------------------
# Specify other directories 
# ------------------------------------------------------------------------

# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS

TAC_ARG_ENABLE_FEATURE(triutils, [Build triutils], TRIUTILS, yes)
AM_CONDITIONAL(BUILD_TRIUTILS, test "X$ac_cv_use_triutils" != "Xno")

# y12m is needed only when AZ_lu is enabled
if test "X$ac_cv_use_aztecoo_azlu" != "Xno"; then
  TAC_ARG_ENABLE_FEATURE(referencey12m, [Build y12m lib included in the distribution], Y12M, no)
  AM_CONDITIONAL(BUILD_Y12M, test "X$ac_cv_use_referencey12m" != "Xno")
# y12m library needed
  if test "X$ac_cv_use_referencey12m" = "Xno"; then
    # Find a valid y12m
    AC_LANG_SAVE
    AC_LANG_FORTRAN77
    AC_SEARCH_LIBS(y12mae,[y12m],,
    AC_MSG_ERROR([[Cannot find y12m library, please specify additional locations using --with-libdirs or type --enable-referencey12m on the configure line in the Trilinos level directory to use the y12m library included with the distribution]]))
    AC_LANG_RESTORE
  fi
fi
# Checks for library functions.
AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
#Removed the following line because it caused problems on janus.
#Also, other libraries check only one function, a standard math
#lib should have sqrt provided that it has pow.
#AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor pow sqrt strchr uname])

# Check for user-specified names for LAPACK and BLAS
TAC_ARG_WITH_LAPACKLIB
TAC_ARG_WITH_BLASLIB

# Define F77_FUNC that will be used to link with Fortran subroutines, e.g., BLAS/LAPACK
AC_F77_WRAPPERS

# find or choose to supply valid blas and lapack libraries
#TAC_ARG_ENABLE_FEATURE(referenceblasandlapack, [Build blas and lapack libs included in the distribution], BLASANDLAPACK, no)
#AM_CONDITIONAL(BUILD_BLAS, test "X$ac_cv_use_referenceblasandlapack" != "Xno")
#if test "X$ac_cv_use_referenceblasandlapack" = "Xno"; then
#  ACX_BLAS([], [AC_MSG_ERROR([[Cannot find blas library, please specify additional locations using --with-libdirs, additional lib names using --with-blaslib or type --with-referenceblasandlapack on the configure line in the Trilinos level directory to use the blas and lapack libraries included in the distribution]])])
  ACX_BLAS([], [AC_MSG_ERROR([[Cannot find blas library, please specify additional locations using --with-libdirs or additional lib names using --with-blaslib]])])
#fi
#TAC_ARG_ENABLE_FEATURE(referencelapack, [Build lapack lib included in the distribution], LAPACK, no)
#AM_CONDITIONAL(BUILD_LAPACK, test "X$ac_cv_use_referencelapack" != "Xno" || test "X$ac_cv_use_referenceblasandlapack" != "Xno")
#if test "X$ac_cv_use_referencelapack" = "Xno" && test "X$ac_cv_use_referenceblasandlapack" = "Xno"; then
#  ACX_LAPACK([], [AC_MSG_ERROR([[Cannot find lapack library, please specify additional locations using --with-libdirs, additional lib names using --with-lapacklib or type --with-referencelapack on the configure line in the Trilinos level directory to use the lapack library included in the distribution]])])
  ACX_LAPACK([], [AC_MSG_ERROR([[Cannot find lapack library, please specify additional locations using --with-libdirs or additional lib names using --with-lapacklib]])])
#fi

# Determine libraries needed for linking with Fortran
AC_F77_LIBRARY_LDFLAGS

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
# Perform substitutions in output files
# ------------------------------------------------------------------------

AC_SUBST(ac_aux_dir)

#-------------------------------------------------------------------------
# Add any flags that are required for external libraries that are to be
# built during the build process.
#-------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Output files
# ------------------------------------------------------------------------
# should have this dependent on if AztecOO is built, or 
# just core Aztec - makefiles in example/
# even more critical

#
#                example/IfpackAztecOO/Makefile
#                example/azoo_iterate_hb/Makefile
#                example/aztec_app/Makefile
#                example/aztec_hb/Makefile
#
AC_CONFIG_FILES([
		Makefile
                src/Makefile
		example/Makefile
                example/simple_epetra_aztecoo/Makefile
                example/AztecOO_UserOpUserMat/Makefile
                example/AztecOO/Makefile
		test/Makefile
		test/AztecOO_UserOpUserMat/Makefile
		test/scripts/Makefile
		test/scripts/daily/Makefile
		test/scripts/daily/serial/Makefile])
AC_OUTPUT

# Bye World!
echo "----------------------------------------"
echo "Finished Running AztecOO Configure Script"
echo "----------------------------------------"
