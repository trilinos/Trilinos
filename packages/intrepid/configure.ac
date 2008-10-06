# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------
# @HEADER
# ************************************************************************
#
#                           Intrepid Package
#                 Copyright (2007) Sandia Corporation
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
# Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
#                    Denis Ridzal (dridzal@sandia.gov).
#
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# Initialization
# ------------------------------------------------------------------------

AC_INIT([intrepid], [1.0], [pbboche@sandia.gov,dridzal@sandia.gov])

echo "----------------------------------------"
echo "Running Intrepid Configure Script"
echo "----------------------------------------"

# ------------------------------------------------------------------------
# This is to protect against accidentally specifying the wrong
# directory with --srcdir.  Any file in that directory will do,
# preferably one that is unlikely to be removed or renamed.
# ------------------------------------------------------------------------
#
AC_CONFIG_SRCDIR([src/local/Intrepid_MultiCell.hpp])



# ------------------------------------------------------------------------
# This a subdirectory containing useful m4 files. The m4 files perform
# various common checks and set common variable values to minimize
# redundant configure work.
# ------------------------------------------------------------------------
#
AC_CONFIG_AUX_DIR(config)



# ------------------------------------------------------------------------
# Configure should create src/Intrepid_config.h from src/Intrepid_config.h.in
# ------------------------------------------------------------------------
#
AM_CONFIG_HEADER(src/Intrepid_config.h:src/Intrepid_config.h.in)


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


# ------------------------------------------------------------------------
# AM_MAINTAINER_MODE turns off maintainer-only makefile targets by
# default, and changes configure to understand a
# --enable-maintainer-mode option. --enable-maintainer-mode turns the
# maintainer-only targets back on. The maintainer-only makefile
# targets permit end users to clean automatically-generated files such
# as configure, which means they have to have autoconf and automake
# installed to repair the damage. AM_MAINTAINER_MODE makes it a bit
# harder for users to shoot themselves in the foot.
# ------------------------------------------------------------------------
#
AM_MAINTAINER_MODE



# ------------------------------------------------------------------------
# Define $build, $host, $target, etc.
# ------------------------------------------------------------------------
#
AC_CANONICAL_TARGET



# ------------------------------------------------------------------------
# Use automake
# ------------------------------------------------------------------------
#  - Required version of automake.
#
AM_INIT_AUTOMAKE(1.10 no-define tar-ustar)

#  - Required version of autoconf.
#
AC_PREREQ(2.61)



# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------
#
TAC_ARG_CONFIG_MPI



# ------------------------------------------------------------------------
# Check for compilers and ranlibs.
# ------------------------------------------------------------------------
#
AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ c++ cxx)
# no need to check for fortran AC_PROG_F77(f77 g77 gfortran f90 xlf90 f95)
AC_PROG_RANLIB



# ------------------------------------------------------------------------
# Check if --with-flags present, prepend any specs to FLAGS
# Check for --with-libs.
# ------------------------------------------------------------------------
#
TAC_ARG_WITH_FLAGS(ccflags, CCFLAGS)
TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
TAC_ARG_WITH_FLAGS(cflags, CFLAGS)
#TAC_ARG_WITH_FLAGS(fflags, FFLAGS)
TAC_ARG_WITH_LIBS
TAC_ARG_WITH_FLAGS(ldflags, LDFLAGS)



# ------------------------------------------------------------------------
# Alternate archiver
# ------------------------------------------------------------------------
#
TAC_ARG_WITH_AR


# ------------------------------------------------------------------------
# MPI link check
# ------------------------------------------------------------------------
TAC_ARG_CHECK_MPI


# ------------------------------------------------------------------------
# Default compiler for header checks.
# Not needed if building within Trilinos.
# ------------------------------------------------------------------------
#
# AC_LANG_CPLUSPLUS


# ------------------------------------------------------------------------
# Checks for Makefile.export related systems
# ------------------------------------------------------------------------
TAC_ARG_ENABLE_EXPORT_MAKEFILES(yes)


# ------------------------------------------------------------------------
# Checks for header files.
# ------------------------------------------------------------------------
#
AC_CHECK_HEADERS([cstdlib stdlib.h], [break],)
AC_CHECK_HEADERS([cstdio stdio.h], [break],)
AC_CHECK_HEADERS([cassert assert.h], [break],)
AC_CHECK_HEADERS([cmath math.h], [break],)
AC_CHECK_HEADERS([iostream iostream.h], [break],)
AC_CHECK_HEADERS([sstream sstream.h], [break],)
AC_CHECK_HEADERS([string string.h], [break],)
AC_CHECK_HEADERS([vector vector.h], [break],)
AC_CHECK_HEADERS([map map.h], [break],)
AC_CHECK_HEADERS([iomanip iomanip.h], [break],)
AC_CHECK_HEADERS([time time.h], [break],)

if test "X$ac_cv_header_cstdlib" != "Xyes" && test "X$ac_cv_header_stdlib_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdlib or stdlib.h.  Intrepid requires at least one or the other.]])
fi

if test "X$ac_cv_header_cstdio" != "Xyes" && test "X$ac_cv_header_stdio_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdio or stdio.h.  Intrepid requires at least one or the other.]])
fi

if test "X$ac_cv_header_cassert" != "Xyes" && test "X$ac_cv_header_assert_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cassert or assert.h.  Intrepid requires at least one or the other.]])
fi

if test "X$ac_cv_header_string" != "Xyes" && test "X$ac_cv_header_string_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find string or string.h.  Intrepid requires at least one or the other.]])
fi

if test "X$ac_cv_header_iostream" != "Xyes" && test "X$ac_cv_header_iostream_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find iostream or iostream.h.  Intrepid requires at least one or the other.]])
fi

if test "X$ac_cv_header_sstream" != "Xyes" && test "X$ac_cv_header_sstream_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find sstream or sstream.h.  Intrepid requires at least one or the other.]])
fi

if test "X$ac_cv_header_cmath" != "Xyes" && test "X$ac_cv_header_math_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cmath or math.h.  Intrepid requires at least one or the other.]])
fi

if test "X$ac_cv_header_vector" != "Xyes" && test "X$ac_cv_header_vector_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find vector or vector.h.  Intrepid requires at least one or the other.]])
fi

if test "X$ac_cv_header_map" != "Xyes" && test "X$ac_cv_header_map_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find map or map.h.  Intrepid requires at least one or the other.]])
fi

if test "X$ac_cv_header_iomanip" != "Xyes" && test "X$ac_cv_header_iomanip_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find iomanip or iomanip.h.  Intrepid requires at least one or the other.]])
fi


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
# External (3rd party) dependencies.
# ------------------------------------------------------------------------
#TAC_ARG_ENABLE_FEATURE_SUB( intrepid, moab, [Enable MOAB interface support.],
#                INTREPID_MOAB, no)
#AM_CONDITIONAL(BUILD_INTREPID_MOAB, [test "X$ac_cv_use_intrepid_moab" != "Xno"])

#if test "X$ac_cv_use_intrepid_moab" != "Xno"; then
#  TODO:  Change to MOAB.
#
#  Check for ARPREC headers.
#
#     AC_MSG_CHECKING(for ARPREC headers (mp/mpreal.h) )
#     AC_TRY_CPP([#include "mp/mpreal.h"],
#        [AC_MSG_RESULT(yes)],
#        [
#         AC_MSG_RESULT(no)
#         echo "-----"
#         echo "Cannot find ARPREC headers."
#         echo "Use --with-incdirs to specify ARPREC include directory."
#         echo "-----"
#         AC_MSG_ERROR(ARPREC headers cannot be found)
#        ])
#
#  Check for ARPREC library.
#
#     AC_MSG_CHECKING(for ARPREC library (libmp.a) )
#     AC_LANG_SAVE
#     AC_LANG_CPLUSPLUS
#     AC_LINK_IFELSE([ AC_LANG_PROGRAM([[#include "mp/mpreal.h"]],
#        [[mp::mp_init(200);]]) ],
#        [AC_MSG_RESULT(yes)],
#        [AC_MSG_RESULT(no)
#                 echo "-----"
#                 echo "Cannot find ARPREC library."
#                 echo "Use --with-ldflags to specify ARPREC library directory."
#                 echo "Use --with-libs to specify ARPREC library (-lmp)."
#                 echo "-----"
#                 AC_MSG_ERROR(ARPREC test cannot link)])
#     AC_LANG_RESTORE
#fi

# ------------------------------------------------------------------------
# Checks if tests and examples should be built
# ------------------------------------------------------------------------

#  #np# - These options can disable the tests and examples of a package.
#  #np# - Packages that do not have tests or examples should #-out the
#  #np# - option(s) that does (do) not apply.

TAC_ARG_ENABLE_FEATURE(tests, [Make tests for all Trilinos packages buildable with 'make tests'], TESTS, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( intrepid, tests, [Make Intrepid tests buildable with 'make tests'], INTREPID_TESTS)
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_intrepid_tests" != "Xno")

TAC_ARG_ENABLE_FEATURE(examples, [Make examples for all Trilinos packages buildable with 'make examples'], EXAMPLES, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( intrepid, examples, [Make Intrepid examples buildable with 'make examples'], NEW_PACKAGE_EXAMPLES)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_intrepid_examples" != "Xno")

# We now build tests and examples through separate make targets, rather than
# during "make".  We still need to conditionally include the test and example
# in SUBDIRS, even though SUB_TEST and SUB_EXAMPLE will never be
# defined, so that the tests and examples are included in the distribution
# tarball.
AM_CONDITIONAL(SUB_TEST, test "X$ac_cv_use_sub_test" = "Xyes")
AM_CONDITIONAL(SUB_EXAMPLE, test "X$ac_cv_use_sub_example" = "Xyes")


# ------------------------------------------------------------------------
# Check to see if the intrepid debug mode has been enabled. In this mode,
# the basic 1-3D linear algebra constructors and operations are
# boundschecked. Other low-level checks are also enabled.
# ------------------------------------------------------------------------
TAC_ARG_ENABLE_FEATURE_SUB( intrepid, debug, [Enable intrepid debug mode],
                INTREPID_DEBUG, no)

# ------------------------------------------------------------------------
# Check to see if FIAT very high order bases have been enabled.
# ------------------------------------------------------------------------
TAC_ARG_ENABLE_FEATURE_SUB( intrepid, fiatvho, [Enable intrepid FIAT very high order bases],
                INTREPID_FIAT_VERY_HIGH_ORDER, no)

# ------------------------------------------------------------------------
# Trilinos Package Dependencies
# Required: Teuchos
# Optional: Epetra, EpetraExt, Sacado, Amesos (w/ Umfpack)
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_CAN_USE_PACKAGE(intrepid, sacado, INTREPID_SACADO, BUILD_SACADO, yes, [Build examples which use Sacado.], [Build examples which use Sacado. Can be overridden with --disable-intrepid-sacado.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(intrepid, epetra, INTREPID_EPETRA, BUILD_EPETRA, yes, [Build examples which use Epetra.], [Build examples which use Epetra. Can be overridden with --disable-intrepid-epetra.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(intrepid, epetraext, INTREPID_EPETRAEXT, BUILD_EPETRAEXT, yes, [Build examples which use EpetraExt.], [Build examples which use EpetraExt. Can be overridden with --disable-intrepid-epetraext.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(intrepid, amesos, INTREPID_AMESOS, BUILD_AMESOS, yes, [Build examples which use Amesos.], [Build examples which use Amesos. Can be overridden with --disable-intrepid-amesos.])

TAC_ARG_INFORM_SUBPACKAGE_USED(amesos, umfpack, INTREPID_AMESOS_UMFPACK, BUILD_AMESOS_UMFPACK, yes, [Build examples which use Amesos interface to Umfpack. Can be overridden with --disable-intrepid-amesos-umfpack.], [Build examples which use Amesos interface to Umfpack.])

# ------------------------------------------------------------------------
# Specify other directories
# ------------------------------------------------------------------------
#
# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS



# ------------------------------------------------------------------------
# Perform substitutions in output files
# ------------------------------------------------------------------------

AC_SUBST(ac_aux_dir)



# ------------------------------------------------------------------------
# List all necessary Makefiles here.
# ------------------------------------------------------------------------
#
AC_CONFIG_FILES([
  Makefile
  Makefile.export.intrepid
  src/Makefile
  test/Makefile
  example/Makefile
  test/Basis/Makefile
  test/Cubature/Makefile
  test/MultiCell/Makefile
  test/RealSpace/Makefile
  test/FieldContainer/Makefile
  test/Basis/F0_QUAD_I1_FEM_DEFAULT/Makefile
  test/Basis/F0_TRI_C1_FEM_DEFAULT/Makefile
  test/Basis/F0_TRI_C2_FEM_DEFAULT/Makefile
  example/Cubature/Makefile
  example/MultiCell/Makefile
  example/RealSpace/Makefile
  example/FieldContainer/Makefile
  example/Basis/Makefile
  example/SparseConv/Makefile
  test/Basis/F0_TRI_C1_FEM_FIAT/Makefile
  test/Basis/F0_TRI_C2_FEM_FIAT/Makefile
  test/Basis/F0_TRI_C3_FEM_FIAT/Makefile
  test/Basis/F0_TRI_C4_FEM_FIAT/Makefile
  test/Basis/F0_TRI_C5_FEM_FIAT/Makefile
  test/Basis/F0_TRI_C6_FEM_FIAT/Makefile
  test/Basis/F0_TRI_C7_FEM_FIAT/Makefile
  test/Basis/F0_TRI_C8_FEM_FIAT/Makefile
  test/Basis/F0_TRI_C9_FEM_FIAT/Makefile
  test/Basis/F0_TRI_C10_FEM_FIAT/Makefile
  test/Basis/F0_TET_C1_FEM_FIAT/Makefile
  test/Basis/F0_TET_C2_FEM_FIAT/Makefile
  test/Basis/F0_TET_C3_FEM_FIAT/Makefile
  test/Basis/F0_TET_C4_FEM_FIAT/Makefile
  test/Basis/F0_TET_C5_FEM_FIAT/Makefile
  test/Basis/F0_TET_C6_FEM_FIAT/Makefile
  test/Basis/F0_TET_C7_FEM_FIAT/Makefile
  test/Basis/F0_TET_C8_FEM_FIAT/Makefile
  test/Basis/F0_TET_C9_FEM_FIAT/Makefile
  test/Basis/F0_TET_C10_FEM_FIAT/Makefile
  test/Basis/F2_TRI_I1_FEM_FIAT/Makefile
  test/Basis/F2_TRI_I2_FEM_FIAT/Makefile
  test/Basis/F2_TRI_I3_FEM_FIAT/Makefile
  test/Basis/F2_TRI_I4_FEM_FIAT/Makefile
  test/Basis/F2_TRI_I5_FEM_FIAT/Makefile
  test/Basis/F2_TRI_I6_FEM_FIAT/Makefile
  test/Basis/F2_TRI_I7_FEM_FIAT/Makefile
  test/Basis/F2_TRI_I8_FEM_FIAT/Makefile
  test/Basis/F2_TRI_I9_FEM_FIAT/Makefile
  test/Basis/F2_TRI_I10_FEM_FIAT/Makefile
  test/Basis/F2_TET_I1_FEM_FIAT/Makefile
  test/Basis/F2_TET_I2_FEM_FIAT/Makefile
  test/Basis/F2_TET_I3_FEM_FIAT/Makefile
  test/Basis/F2_TET_I4_FEM_FIAT/Makefile
  test/Basis/F2_TET_I5_FEM_FIAT/Makefile
  test/Basis/F2_TET_I6_FEM_FIAT/Makefile
  test/Basis/F2_TET_I7_FEM_FIAT/Makefile
  test/Basis/F2_TET_I8_FEM_FIAT/Makefile
  test/Basis/F2_TET_I9_FEM_FIAT/Makefile
  test/Basis/F2_TET_I10_FEM_FIAT/Makefile
  test/Basis/F1_TRI_I1_FEM_FIAT/Makefile
  test/Basis/F1_TRI_I2_FEM_FIAT/Makefile
  test/Basis/F1_TRI_I3_FEM_FIAT/Makefile
  test/Basis/F1_TRI_I4_FEM_FIAT/Makefile
  test/Basis/F1_TRI_I5_FEM_FIAT/Makefile
  test/Basis/F1_TRI_I6_FEM_FIAT/Makefile
  test/Basis/F1_TRI_I7_FEM_FIAT/Makefile
  test/Basis/F1_TRI_I8_FEM_FIAT/Makefile
  test/Basis/F1_TRI_I9_FEM_FIAT/Makefile
  test/Basis/F1_TRI_I10_FEM_FIAT/Makefile
  test/Basis/F1_TET_I1_FEM_FIAT/Makefile
  test/Basis/F1_TET_I2_FEM_FIAT/Makefile
  test/Basis/F1_TET_I3_FEM_FIAT/Makefile
  test/Basis/F1_TET_I4_FEM_FIAT/Makefile
  test/Basis/F1_TET_I5_FEM_FIAT/Makefile
  test/Basis/F1_TET_I6_FEM_FIAT/Makefile
  test/Basis/F1_TET_I7_FEM_FIAT/Makefile
  test/Basis/F1_TET_I8_FEM_FIAT/Makefile
  test/Basis/F1_TET_I9_FEM_FIAT/Makefile
  test/Basis/F1_TET_I10_FEM_FIAT/Makefile
  test/Basis/FIAT/Makefile
])

AC_OUTPUT()

echo "------------------------------------------"
echo "Finished Running Intrepid Configure Script"
echo "------------------------------------------"
