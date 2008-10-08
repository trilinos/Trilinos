# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
# 
#        Phalanx: A Partial Differential Equation Field Evaluation 
#       Kernel for Flexible Management of Complex Dependency Chains
#                  Copyright (2008) Sandia Corporation
# 
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
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
# 
# Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
# National Laboratories.
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# Initialization 
# ------------------------------------------------------------------------



# This must be the first line in configure.ac.
# Optional 3rd argument is email address for bugs.

#  #np# - package name, version number, and e-mail address below
AC_INIT(phalanx, 1.0, rppawlo@sandia.gov)

# Hello World!
#  #auto np# - package name below
echo "----------------------------------------"
echo "Running Phalanx Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.  Any file in that directory will do,
# preferably one that is unlikely to be removed or renamed.

# #auto np# - file name below
AC_CONFIG_SRCDIR([src/Phalanx_FieldTag.hpp])

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

#
# #np# This a a subdirectory of Trilinos/packages/new_package that contains
# #np# a number of useful m4 files.  Any package that uses these m4 files
# #np# will have to add them to the repository.  The m4 files perform various
# #np# common checks and set common variable values to minimize redundant work.
#
AC_CONFIG_AUX_DIR(config)
#  #auto np# - Change file names in next line 
# Configure should create src/New_Package_config.h from src/New_Package_config.h.in

AM_CONFIG_HEADER(src/Phalanx_config.hpp:src/Phalanx_config.hpp.in)

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

#  - Required version of automake.
AM_INIT_AUTOMAKE(1.10 no-define tar-ustar)

# Specify required version of autoconf.

AC_PREREQ(2.61)

# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI

#  #np# - can eliminate compiler checks below if your package does not use the
#         language corresponding to the check.  Please note that if you use
#	  F77_FUNC to determine Fortran name mangling, you should not remove
#	  the Fortran compiler check or the check for Fortran flags.  Doing
#	  so will prevent the detection of the proper name mangling in some
#	  cases.
# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ c++ cxx)
##AC_PROG_F77(f77 g77 gfortran f90 xlf90 f95)
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Checks for extra flags
# #np# - You can pull flag checks for languages not used in package.  Please 
#	  note that if you use F77_FUNC to determine Fortran name mangling, 
#	  you should not remove the Fortran compiler check or the check for 
#	  Fortran flags.  Doing so will prevent the detection of the proper 
#	  name mangling in some cases.
# ------------------------------------------------------------------------
# Check if --with-flags present, prepend any specs to FLAGS

TAC_ARG_WITH_FLAGS(ccflags, CCFLAGS)
TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
TAC_ARG_WITH_FLAGS(cflags, CFLAGS)
##TAC_ARG_WITH_FLAGS(fflags, FFLAGS)
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

#  #np# - Customize the checks for headers needed in package being added.
#NOTE: This section still needs to be updated for new_package, it is still 
# code from Epetra
# ------------------------------------------------------------------------
# Checks for header files - used in Epetra_ConfigDefs.h 
# ------------------------------------------------------------------------

## Phalanx requires:
## algorithm
## cmath
## cstddef
## iomanip
## iostream
## iterator
## map
## sstream
## string
## typeinfo
## vector

AC_CHECK_HEADERS([algorithm algo.h algorithm.h], [break],)
AC_CHECK_HEADERS([cmath math.h], [break],)
AC_CHECK_HEADERS([cstddef stddef.h], [break],)
AC_CHECK_HEADERS([fstream fstream.h], [break],)
AC_CHECK_HEADERS([iomanip iomanip.h], [break],)
AC_CHECK_HEADERS([iostream iostream.h], [break],)
AC_CHECK_HEADERS([iterator iterator.h], [break],)
AC_CHECK_HEADERS([map map.h], [break],)
AC_CHECK_HEADERS([sstream sstream.h], [break],)
AC_CHECK_HEADERS([string string.h], [break],)
AC_CHECK_HEADERS([typeinfo typeinfo.h], [break],)
AC_CHECK_HEADERS([vector vector.h], [break],)

if test "X$ac_cv_header_algorithm" != "Xyes" && test "X$ac_cv_header_algo_h" != "Xyes" && test "X$ac_cv_header_algorithm_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find algorithm, algo.h, or algorithm.h.  Phalanx requires at least one of these.]])
fi

if test "X$ac_cv_header_cmath" != "Xyes" && test "X$ac_cv_header_math_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cmath or math.h.  Phalanx requires at least one or the other.]])
fi

if test "X$ac_cv_header_cstddef" != "Xyes" && test "X$ac_cv_header_stddef_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstddef or stddef.h.  Phalanx requires at least one or the other.]])
fi

if test "X$ac_cv_header_fstream" != "Xyes" && test "X$ac_cv_header_fstream_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find fstream or fstream.h.  Phalanx requires at least one or the other.]])
fi

if test "X$ac_cv_header_iomanip" != "Xyes" && test "X$ac_cv_header_iomanip_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find iomanip or iomanip.h.  Phalanx requires at least one or the other.]])
fi

if test "X$ac_cv_header_iostream" != "Xyes" && test "X$ac_cv_header_iostream_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find iostream or iostream.h.  Phalanx requires at least one or the other.]])
fi

if test "X$ac_cv_header_iterator" != "Xyes" && test "X$ac_cv_header_iterator_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find iterator or iterator.h.  Phalanx requires at least one or the other.]])
fi

if test "X$ac_cv_header_map" != "Xyes" && test "X$ac_cv_header_map_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find map or map.h.  Phalanx requires at least one or the other.]])
fi

if test "X$ac_cv_header_sstream" != "Xyes" && test "X$ac_cv_header_sstream_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find sstream or sstream.h.  Phalanx requires at least one or the other.]])
fi

if test "X$ac_cv_header_string" != "Xyes" && test "X$ac_cv_header_string_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find string or string.h.  Phalanx requires at least one or the other.]])
fi

if test "X$ac_cv_header_typeinfo" != "Xyes" && test "X$ac_cv_header_typeinfo_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find typeinfo or typeinfo.h.  Phalanx requires at least one or the other.]])
fi

if test "X$ac_cv_header_vector" != "Xyes" && test "X$ac_cv_header_vector_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find vector or vector.h.  Phalanx requires at least one or the other.]])
fi

# ------------------------------------------------------------------------
# Checks for Makefile.export related systems
# ------------------------------------------------------------------------
TAC_ARG_ENABLE_EXPORT_MAKEFILES(yes)

# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------

## Debug Mode
AC_ARG_ENABLE(debug,
[AC_HELP_STRING([--enable-debug],[Trilinos debug mode flag.])],
[TRILINOS_DEBUG_MODE=$enableval],
[TRILINOS_DEBUG_MODE=no]
)

AC_ARG_ENABLE(phalanx-debug,
[AC_HELP_STRING([--enable-phalanx-debug],[Overrides Trilinos debug mode for phalanx package only.])],
[PHALANX_DEBUG_MODE=$enableval],
[PHALANX_DEBUG_MODE=TRILINOS_DEBUG_MODE]
)

AC_MSG_CHECKING(whether phalanx is being built in debug mode)
if test "X$PHALANX_DEBUG_MODE" != "Xno"; then
  AC_DEFINE([HAVE_PHALANX_DEBUG],1,[Define to build Phalanx in debug mode.])
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
fi

# ------------------------------------------------------------------------
# Checks if tests and examples should be built
# ------------------------------------------------------------------------

#  #np# - These options can disable the tests and examples of a package.
#  #np# - Packages that do not have tests or examples should #-out the 
#  #np# - option(s) that does (do) not apply.

TAC_ARG_ENABLE_FEATURE(tests, [Make tests for all Trilinos packages buildable with 'make tests'], TESTS, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( phalanx, tests, [Make Phalanx tests buildable with 'make tests'], PHALANX_TESTS)
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_phalanx_tests" != "Xno")

TAC_ARG_ENABLE_FEATURE(examples, [Make examples for all Trilinos packages buildable with 'make examples'], EXAMPLES, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( phalanx, examples, [Make Phalanx examples buildable with 'make examples'], PHALANX_EXAMPLES)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_phalanx_examples" != "Xno")

#We now build tests and examples through separate make targets, rather than
#during "make".  We still need to conditionally include the test and example
#in SUBDIRS, even though SUB_TEST and SUB_EXAMPLE will never be
#defined, so that the tests and examples are included in the distribution
#tarball.
AM_CONDITIONAL(SUB_TEST, test "X$ac_cv_use_sub_test" = "Xyes")
AM_CONDITIONAL(SUB_EXAMPLE, test "X$ac_cv_use_sub_example" = "Xyes")

TAC_ARG_ENABLE_FEATURE(libcheck, [Check for some third-party libraries.  (Cannot be disabled unless tests and examples are also disabled.)], LIBCHECK, yes)

# ------------------------------------------------------------------------
# Specify other directories 
# ------------------------------------------------------------------------

# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS

# #np# - Yet another opportunity to remove code if you aren't
# using Fortran
# Define F77_FUNC that will be used to link with Fortran subroutines. - trash WORKGXX 
##AC_F77_WRAPPERS

# ------------------------------------------------------------------------
# Checks for required and optional libraries
# ------------------------------------------------------------------------

## Require the Trilinos Libraries: Teuchos and Sacado
## No checks needed - these will automatically be turned on at the top level of Trilinos.  But just to make sure...
TRILINOS_HEADER(teuchos,Teuchos_config.h,src)
TRILINOS_HEADER(teuchos,Teuchos_Flops.hpp,src)
TRILINOS_HEADER(sacado,Sacado_config.h,src)
TRILINOS_HEADER(sacado,Sacado_mpl_at.hpp,src/mpl)

## Require the Boost MPL Library headers only (don't need a check 
## for the actual library).

# The command below fails to pull in --with-incdirs flags.  
##AC_CHECK_HEADERS([boost/mpl/at.h],,AC_MSG_ERROR([[Cannot find the Boost MPL library.  Phalanx requires Boost be installed on your machine.  Use --with-incdirs=-I<DIR> to specify the install directory for Boost.]]))

## Optionally require the TVMET Library headers only (don't need a check 
## for the actual library).
AC_ARG_WITH(tvmet,
[AC_HELP_STRING([--with-tvmet],[Enable TVMET support])],
[
  if test "X${withval}" != "Xyes"; then
     MF_PATH=${withval}
  fi
  HAVE_TVMET=yes
],
[HAVE_TVMET=no]
)
AC_MSG_CHECKING(whether tvmet is enabled)
AC_MSG_RESULT([${HAVE_TVMET}])
if test "X${HAVE_TVMET}" != "Xno"; then

  AC_PREPROC_IFELSE([AC_LANG_SOURCE(
  [[
#undef PACKAGE
#undef PACKAGE_NAME
#undef PACKAGE_BUGREPORT
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION
#include "tvmet/Vector.h"
  ]])],
    [ have_tvmet_header=yes
    ],[have_tvmet_header=no])

  AC_MSG_CHECKING(for tvmet header Vector.h)
  AC_MSG_RESULT([${have_tvmet_header}])
  if test $have_tvmet_header = no; then
    AC_MSG_ERROR([Cannot find tvmet header tvmet/Vector.h.  Please add path to tvmet include directory using --with-incdirs flag.])
  fi

  AC_DEFINE([HAVE_PHALANX_TVMET],1,[Define if tvmet is available])
fi
AM_CONDITIONAL(BUILD_TVMET, test X${BUILD_TVMET} = Xyes)

## Optional dependence on Epetra for example problem
TAC_ARG_ENABLE_CAN_USE_PACKAGE(phalanx, epetra, PHALANX_EPETRA, BUILD_PHALANX_EPETRA, yes, [Build epetra support for phalanx examples.], [Implicitly build epetra support for phalanx examples. Can be overridden with --<en|dis>able-phalanx-epetra.])
## Optional dependence on Ifpack for example problem
TAC_ARG_ENABLE_CAN_USE_PACKAGE(phalanx, ifpack, PHALANX_IFPACK, BUILD_PHALANX_IFPACK, yes, [Build ifpack support for phalanx examples.], [Implicitly build ifpack support for phalanx examples. Can be overridden with --<en|dis>able-phalanx-ifpack.])
## Optional dependence on Belos for example problem
TAC_ARG_ENABLE_CAN_USE_PACKAGE(phalanx, belos, PHALANX_BELOS, BUILD_PHALANX_BELOS, yes, [Build belos support for phalanx examples.], [Implicitly build belos support for phalanx examples. Can be overridden with --<en|dis>able-phalanx-belos.])


# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

# If tests, examples and libcheck are disabled, we don't have to check
# for these libraries.

# #np# -
# If a package does not have tests or examples, the corresponding check(s)
# should be pulled out of the "if" statement below.
#if test "X$ac_cv_use_new_package_examples" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then
if test "X$ac_cv_use_phalanx_tests" != "Xno" || test "X$ac_cv_use_phalanx_examples" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then

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
# Checks for linker characteristics
# ------------------------------------------------------------------------

# Determine libraries needed for linking with Fortran
AC_F77_LIBRARY_LDFLAGS


#  ----------------------------------------------------
#  Needs this because Dec handles templates differently
#  ----------------------------------------------------
#  #np# - Only uncomment the following lines if the package
#         creates a subdirectory called cxx_repository for
#         template object files on DEC.  (HINT: try to build on 
#         DEC/CPLANT and see if template-related errors occur)

# -- DEC CXX --
#if test "X$TRILINOS_TEMPLATE_OBJS" = "X"; then
#  case $target in
#    alpha*)
#      if test $CXX = cxx; then
#        TRILINOS_TEMPLATE_OBJS="cxx_repository/*.o"
#      fi
#    ;;
#  esac
#fi
                                                                  
#AC_MSG_CHECKING([for template-specific object files to be linked])if test "X$TRILINOS_TEMPLATE_OBJS" = "X"; then
# AC_MSG_RESULT(none)
#else
# AC_MSG_RESULT($TRILINOS_TEMPLATE_OBJS)
#fi
                                                                  
#AC_SUBST(TRILINOS_TEMPLATE_OBJS)

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
		Makefile.export.phalanx
		doc/Makefile
		src/Makefile
		example/Makefile
		example/EnergyFlux/Makefile
		example/MultiDimensionalArray/Makefile
		example/FEM_Nonlinear/Makefile
		example/FEM_Nonlinear_VBR/Makefile
		test/Makefile
		test/Utilities/Makefile
		test/Allocators/Makefile
		test/DataContainer/Makefile
		test/DataLayout/Makefile
		test/EvaluatorManager/Makefile
		test/FieldTag/Makefile
		test/Field/Makefile
		test/Performance/Makefile
		test/Performance/AlgebraicTypes/Makefile
		test/Performance/Allocators/Makefile
		test/Performance/BracketOperator/Makefile
		])

AC_OUTPUT()

# Bye World!
#  #auto np# change package name
echo "---------------------------------------------"
echo "Finished Running Phalanx Configure Script"
echo "---------------------------------------------"
