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

#  #np# - package name, version number, and e-mail address below
AC_INIT(anasazi, 2.2, hkthorn@sandia.gov)

# Hello World!
#  #np# - package name below
echo "----------------------------------------"
echo "Running Anasazi Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.  Any file in that directory will do,
# preferably one that is unlikely to be removed or renamed.

# #np# - file name below
AC_CONFIG_SRCDIR([src/AnasaziMultiVec.hpp])

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

#
# #np# This a a subdirectory of Trilinos/packages/new_package that contains
# #np# a number of useful m4 files.  Any package that uses these m4 files
# #np# will have to add them to the repository.  The m4 files perform various
# #np# common checks and set common variable values to minimize redundant work.
#
AC_CONFIG_AUX_DIR(config)
#  #np# - Change file names in next line 
# Configure should create src/New_Package_config.h from src/New_Package_config.h.in

AM_CONFIG_HEADER(src/Anasazi_config.h:src/Anasazi_config.h.in)

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

# Specify required version of automake. 
AM_INIT_AUTOMAKE(1.10 no-define tar-ustar)

# Specify required version of autoconf.

AC_PREREQ(2.61)

# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI

#  #np# - can eliminate compiler checks below if your package does not use the
#         language corresponding to the check
# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ gcc c++ cxx)
#AC_PROG_F77(f77 g77 f90 xlf90)
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Checks for extra flags
# #np# - You can pull flag checks for languages not used in package 
# ------------------------------------------------------------------------
# Check if --with-flags present, prepend any specs to FLAGS

TAC_ARG_WITH_FLAGS(ccflags, CCFLAGS)
TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
TAC_ARG_WITH_FLAGS(cflags, CFLAGS)
#TAC_ARG_WITH_FLAGS(fflags, FFLAGS)
TAC_ARG_WITH_LIBS
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
# Checks for Makefile.export related systems
# ------------------------------------------------------------------------
TAC_ARG_ENABLE_EXPORT_MAKEFILES(yes)

#  #np# - Customize the checks for headers needed in package being added.
#NOTE: This section still needs to be updated for new_package, it is still 
# code from Epetra
# ------------------------------------------------------------------------
# Checks for header files - used in Anasazi_ConfigDefs.h 
# ------------------------------------------------------------------------

AC_CHECK_HEADERS([cstdlib stdlib.h], [break],)
AC_CHECK_HEADERS([cstdio stdio.h], [break],)
# AC_CHECK_HEADERS([cassert assert.h], [break],)
AC_CHECK_HEADERS([cmath math.h], [break],)
AC_CHECK_HEADERS([iomanip iomanip.h], [break],)
AC_CHECK_HEADERS([iostream iostream.h], [break],)
AC_CHECK_HEADERS([iterator iterator.h], [break],)
AC_CHECK_HEADERS([string string.h], [break],)
AC_CHECK_HEADERS([vector vector.h], [break],)
AC_CHECK_HEADERS([complex complex.h], [break],)
AC_CHECK_HEADERS([algorithm algo.h algorithm.h], [break],)
AC_CHECK_HEADERS([cctype cctype.h], [break],)
AC_CHECK_HEADERS([stdexcept stdexcept.h], [break],)
AC_CHECK_HEADERS([numeric algo.h],[break],)

if test "X$ac_cv_header_cstdlib" != "Xyes" && test "X$ac_cv_header_stdlib_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdlib or stdlib.h.  Anasazi requires at least one or the other.]])
fi

if test "X$ac_cv_header_cstdio" != "Xyes" && test "X$ac_cv_header_stdio_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdio or stdio.h.  Anasazi requires at least one or the other.]])
fi

#if test "X$ac_cv_header_cassert" != "Xyes" && test "X$ac_cv_header_assert_h" != "Xyes"; then
#  AC_MSG_ERROR([[Cannot find cassert or assert.h.  Anasazi requires at least one or the other.]])
#fi

if test "X$ac_cv_header_string" != "Xyes" && test "X$ac_cv_header_string_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find string or string.h.  Anasazi requires at least one or the other.]])
fi

if test "X$ac_cv_header_iterator" != "Xyes" && test "X$ac_cv_header_iterator_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find iterator or iterator.h.  Anasazi requires at least one or the other.]])
fi

if test "X$ac_cv_header_vector" != "Xyes" && test "X$ac_cv_header_vector_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find vector or vector.h.  Anasazi requires at least one or the other.]])
fi

if test "X$ac_cv_header_complex" != "Xyes" && test "X$ac_cv_header_complex_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find complex or complex.h.  Anasazi requires at least one or the other.]])
fi

if test "X$ac_cv_header_algorithm" != "Xyes" && test "X$ac_cv_header_algorithm_h" != "Xyes" && test "X$ac_cv_header_algo_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find algorithm or algorithm.h.  Anasazi requires at least one or the other.]])
fi

if test "X$ac_cv_header_cctype" != "Xyes" && test "X$ac_cv_header_cctype_h" != "Xyes" && test "X$ac_cv_header_algo_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cctype or cctype.h.  Anasazi requires at least one or the other.]])
fi

if test "X$ac_cv_header_numeric" != "Xyes" && test "X$ac_cv_header_algo_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find numeric or algo.h.  Anasazi requires at least one or the other.]])
fi

if test "X${host}" != "Xtflop"; then
  if test "X$ac_cv_header_cmath" != "Xyes" && test "X$ac_cv_header_math_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find cmath or math.h.  Anasazi requires at least one or the other.]])
  fi
else
  # TFLOP
  if test "X$ac_cv_header_iomanip" != "Xyes" && test "X$ac_cv_header_iomanip_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find iomanip or iomanip.h.  Anasazi requires at least one or the other.]])
  fi
fi  

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE

# ------------------------------------------------------------------------
# Checks for compiler characteristics - Anasazi_ConfigDefs.h uses these, 
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
# a cpp #define command and one or more source or header files  - ie 
# conditional compilation which requires a change to the makefile.
# If configure is called with "--enable-newp_swahili", it will place 
# a #define HAVE_NEWP_SWAHILI in New_Package_config.h and add source 
# and header files as specified in Makefile.am 
#
# To add an option to your package, you will need to modify this file
# (i.e. configure.ac), and src/Makefile.am.  A bootstrap is also required.
# Look for SWAHILI in both files listed above.
#
# TAC_ARG_ENABLE_OPTION performs several tasks, including:
#   1)  Setting a help string
#   2)  Printing a message to standard out
#   3)  Defining HAVE_NEWP_SWAHILI (For use in New_Package_config.h)
#
# AM_CONDITIONAL makes HAVE_NEWP_SWAHILI available for use in Makefile.am
# ------------------------------------------------------------------------
#
#TAC_ARG_ENABLE_OPTION(newp_swahili, [enable swahili support], NEWP_SWAHILI, no)
#AM_CONDITIONAL(HAVE_NEWP_SWAHILI, [test "X$ac_cv_use_newp_swahili" != "Xno"])

# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------

#
#  #np# - These options are specific to Anasazi and are left only as examples.
#         If uncommented, the following lines would define HAVE_FATAL_MESSAGES
#         in the (package)_config.h file.  Developers should use a 
#         more specific name such as 'new_package-fatal_messages', rather
#         than 'fatal_messages'.
#TAC_ARG_ENABLE_OPTION(fatal_messages, [enable printing of fatal error messages], 
#                   FATAL_MESSAGES, yes)
#TAC_ARG_ENABLE_OPTION(warning_messages, [enable printing of warning error messages], 
#                   WARNING_MESSAGES, no)
#TAC_ARG_ENABLE_OPTION(format_io, [Use package-defined IO formatting for output], 
#                   FORMAT_IO, yes)

#  #np# - These are specific to Anasazi and are left only as examples.
#         If uncommented, the following line would define HAVE_ZOLTAN in
#         the (package)_config.h file.  Developers should use a more
#         specific name such as 'new_package-zoltan', rather than 'zoltan'.
# TAC_ARG_WITH_PACKAGE(zoltan, [Enable Zoltan interface support], ZOLTAN, no)
# TAC_ARG_WITH_PACKAGE(rtop, [Enable RTOp interface support], RTOP, no)

# ------------------------------------------------------------------------
# Checks if examples should be built or which ones should be built
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE(examples, [Make examples for all Trilinos packages buildable with 'make examples'], EXAMPLES, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( anasazi, examples, [Make Anasazi examples buildable with 'make examples'], ANASAZI_EXAMPLES)
AM_CONDITIONAL(BUILD_EXAMPLES, [test "X$ac_cv_use_anasazi_examples" != "Xno"])

TAC_ARG_ENABLE_FEATURE(tests, [Make tests for all Trilinos packages buildable with 'make tests'], TESTS, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( anasazi, tests, [Make Anasazi tests buildable with 'make tests'], ANASAZI_TESTS)
AM_CONDITIONAL(BUILD_TESTS, [test "X$ac_cv_use_anasazi_tests" != "Xno"])

#We now build tests and examples through separate make targets, rather than
#during "make".  We still need to conditionally include the test and example
#in SUBDIRS, even though SUB_TEST and SUB_EXAMPLE will never be
#defined, so that the tests and examples are included in the distribution
#tarball.
AM_CONDITIONAL(SUB_TEST, test "X$ac_cv_use_sub_test" = "Xyes")
AM_CONDITIONAL(SUB_EXAMPLE, test "X$ac_cv_use_sub_example" = "Xyes")

TAC_ARG_ENABLE_FEATURE(libcheck, [Check for some third-party libraries.  (Cannot be disabled unless tests and examples are also disabled.)], LIBCHECK, yes)

# ------------------------------------------------------------------------
# Trilinos Package Dependencies
# Required: Teuchos
# Optional: Epetra, Tpetra, Thyra, Triutils, Belos, AztecOO, and Ifpack
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE_SUB(teuchos, complex, [Build examples/tests which use complex scalar type], TEUCHOS_COMPLEX, yes) 
AM_CONDITIONAL(BUILD_TEUCHOS_COMPLEX, [test "X$ac_cv_use_teuchos_complex" != "Xno"])

TAC_ARG_ENABLE_FEATURE_SUB(anasazi, experimental, [Build examples/tests which use experimental code], ANASAZI_EXPERIMENTAL, no) 
AM_CONDITIONAL(BUILD_ANASAZI_EXPERIMENTAL, [test "X$ac_cv_use_anasazi_experimental" != "Xno"])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(anasazi, triutils, ANASAZI_TRIUTILS, BUILD_TRIUTILS, yes, [Build examples which use Triutils.], [Build examples which use Triutils. Can be overridden with --disable-anasazi-triutils.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(anasazi, epetra, ANASAZI_EPETRA, BUILD_ANASAZI_EPETRA, yes, [Build examples which use Epetra.], [Build examples which use Epetra. Can be overridden with --disable-anasazi-epetra.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(anasazi, epetraext, ANASAZI_EPETRAEXT, BUILD_EPETRAEXT, yes, [Build examples which use EpetraExt.], [Build examples which use EpetraExt. Can be overridden with --disable-anasazi-epetraext.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(anasazi, belos, ANASAZI_BELOS, BUILD_BELOS, yes, [Build examples which use Belos. Can be overridden with --disable-anasazi-belos.], [Build examples which use Belos.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(anasazi, aztecoo, ANASAZI_AZTECOO, BUILD_AZTECOO, yes, [Build examples which use AztecOO. Can be overridden with --disable-anasazi-aztecoo.], [Build examples which use AztecOO.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(anasazi, amesos, ANASAZI_AMESOS, BUILD_AMESOS, yes, [Build examples which use Amesos. Can be overridden with --disable-anasazi-amesos.], [Build examples which use Amesos.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(anasazi, ifpack, ANASAZI_IFPACK, BUILD_IFPACK, yes, [Build examples which use Ifpack. Can be overridden with --disable-anasazi-ifpack.], [Build examples which use Ifpack.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(anasazi, thyra, ANASAZI_THYRA, BUILD_ANASAZI_THYRA, yes, [Build code which uses Thyra. Can be overridden with --disable-anasazi-thyra.], [Build code which uses Thyra.])

TAC_ARG_INFORM_SUBPACKAGE_USED(epetra, thyra, EPETRA_THYRA, BUILD_EPETRA_THYRA, yes, [Build examples which use Epetra-Thyra adapters. Can be overridden with --disable-epetra-thyra.], [Build examples which use Epetra-Thyra adapters.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(anasazi, tpetra, ANASAZI_TPETRA, BUILD_ANASAZI_TPETRA, no, [Build examples which use Tpetra adapters.  Can be overridden with --disable-anasazi-tpetra.], [Build examples which use Tpetra adapters.])

# ------------------------------------------------------------------------
# Specify other directories 
# ------------------------------------------------------------------------

# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS

# #np# - Yet another opportunity to remove code if you aren't
# using Fortran
# Define F77_FUNC that will be used to link with Fortran subroutines.
#AC_F77_WRAPPERS

# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

# If tests, examples and libcheck are disabled, we don't have to check
# for these libraries.
if test "X$ac_cv_use_tests" != "Xno" || test "X$ac_cv_use_examples" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then

# #np# - 
# To add any libraries that you want to check for, replace pow below 
# with the name of an entry point in that library.
# checks for libraries now
dnl Replace `main' with a function in -lm:
#  #np# These lines can be removed if the math lib is not needed.
AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))

fi
# end of the list of libraries that don't need to be checked for if
# tests and examples are disabled.

# ------------------------------------------------------------------------
# Checks for compiler characteristics
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for linker characteristics
# ------------------------------------------------------------------------

# Determine libraries needed for linking with Fortran
#AC_F77_LIBRARY_LDFLAGS


#  ----------------------------------------------------
#  Needs this because Dec handles templates differently
#  ----------------------------------------------------
#  #np# - Only uncomment the following lines if the package
#         creates a subdirectory called cxx_repository for
#         template object files on DEC.  (HINT: try to build on 
#         DEC/CPLANT and see if template-related errors occur)

# -- DEC CXX --
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
#  You will need to change AC_CONFIG_FILES below and Makefile.am
#  to add a new directory.
AC_CONFIG_FILES([
		Makefile 
		Makefile.export.anasazi
		src/Makefile
		test/Makefile
		test/ARPACKExamples/Makefile
		test/MVOPTester/Makefile
 		test/IRTR/Makefile
		test/LOBPCG/Makefile
		test/BlockDavidson/Makefile
		test/BlockKrylovSchur/Makefile
		test/OrthoManager/Makefile
		test/StatusTest/Makefile
		test/SortManager/Makefile
		test/scripts/Makefile
		test/scripts/daily/Makefile
		test/scripts/daily/serial/Makefile
		test/scripts/daily/mpi/Makefile
		example/Makefile
		epetra/Makefile
		epetra/src/Makefile
		epetra/util/Makefile
		epetra/util/ModeLaplace/Makefile
		epetra/example/Makefile
		epetra/example/BlockDavidson/Makefile
		epetra/example/BlockKrylovSchur/Makefile
		epetra/example/LOBPCG/Makefile
		epetra/example/MVOPTester/Makefile
		epetra/test/Makefile
		epetra/test/MVOPTester/Makefile
		epetra/test/ModalSolverUtils/Makefile
		epetra/test/BlockDavidson/Makefile
		epetra/test/BlockKrylovSchur/Makefile
		epetra/test/LOBPCG/Makefile
		epetra/test/IRTR/Makefile
		epetra/test/OrthoManager/Makefile
		thyra/Makefile
		thyra/src/Makefile
		thyra/example/Makefile
		thyra/test/Makefile
		thyra/example/MVOPTester/Makefile
		thyra/test/MVOPTester/Makefile
 		thyra/test/IRTR/Makefile
		thyra/test/LOBPCG/Makefile
		thyra/test/BlockKrylovSchur/Makefile
		thyra/test/BlockDavidson/Makefile
		tpetra/Makefile
		tpetra/src/Makefile
		tpetra/test/Makefile
		tpetra/test/MVOPTester/Makefile
		tpetra/test/BlockKrylovSchur/Makefile
		tpetra/test/IRTR/Makefile
		tpetra/test/LOBPCG/Makefile
		tpetra/test/BlockDavidson/Makefile
		])

AC_OUTPUT()

# Bye World!
#  #np# change package name
echo "----------------------------------------"
echo "Finished Running Anasazi Configure Script"
echo "----------------------------------------"
