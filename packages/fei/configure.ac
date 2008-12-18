# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
#
#              Finite Element Interface to Linear Solvers
#                 Copyright (2004) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# Questions? Contact Alan Williams (william@sandia.gov)
#
# ************************************************************************
# @HEADER

AC_INIT(FEI, 2.23.02, william@sandia.gov)

echo "----------------------------------------"
echo "Running FEI Configure Script"
echo "----------------------------------------"

AC_CONFIG_AUX_DIR(config)

AM_INIT_AUTOMAKE([no-define])

AC_CONFIG_HEADERS([FEI_config.h:FEI_config.h.in])

AM_MAINTAINER_MODE

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
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

AC_PROG_CXX(CC cxx ciCC g++ gcc c++)
AC_PROG_RANLIB

AC_LANG_CPLUSPLUS

AH_TOP([
#ifndef _FEI_config_h_
#define _FEI_config_h_
])

AH_BOTTOM([
#endif
])

# Autoheader always insists on creating #defines for the following
# 5 macros, which can occasionally cause clashes with macros of the same
# name that originate in other auto-conf'd packages. We'll use the
# AH_VERBATIM command to put '#if 0' blocks around all of these names.

AH_VERBATIM([PACKAGE_VERSION], [
/* we use '#if 0' to avoid defining PACKAGE_VERSION */
#if 0
#undef PACKAGE_VERSION
#endif])
AH_VERBATIM([PACKAGE_BUGREPORT], [
/* we use '#if 0' to avoid defining PACKAGE_BUGREPORT */
#if 0
#undef PACKAGE_BUGREPORT
#endif])
AH_VERBATIM([PACKAGE_STRING], [
/* we use '#if 0' to avoid defining PACKAGE_STRING */
#if 0
#undef PACKAGE_STRING
#endif])
AH_VERBATIM([PACKAGE_TARNAME], [
/* we use '#if 0' to avoid defining PACKAGE_TARNAME */
#if 0
#undef PACKAGE_TARNAME
#endif])
AH_VERBATIM([PACKAGE_NAME], [
/* we use '#if 0' to avoid defining PACKAGE_NAME */
#if 0
#undef PACKAGE_NAME
#endif])

TAC_ARG_WITH_FLAGS(cppflags, CPPFLAGS)
TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
TAC_ARG_WITH_LIBS
TAC_ARG_WITH_FLAGS(ldflags, LDFLAGS)


#
# We now check for certain system headers that may not always be found, and
# which we can do without.
#
# We don't bother checking for headers that are always present, such as
# <cstdlib>, <map>, <vector>, etc. If one of those isn't available, then
# FEI simply can't be built.

AC_CHECK_HEADERS([time.h], [],
                 [AC_DEFINE([HAVE_NO_TIME_H], 1, [<time.h> not available])],
                 [#include <stdlib.h>])

AC_CHECK_HEADERS([iosfwd], [],
                 [AC_DEFINE([HAVE_NO_IOSFWD], 1, [<iosfwd> not available])],
                 [#include <cstdlib>])

AC_CHECK_HEADERS([iomanip], [],
                 [AC_DEFINE([HAVE_NO_IOMANIP], 1, [<iomanip> not available])],
                 [#include <cstdlib>])

AC_CHECK_HEADERS([iostream], [],
                 [AC_DEFINE([HAVE_NO_IOSTREAM], 1, [<iostream> not available])],
                 [#include <cstdlib>])

AC_CHECK_HEADERS([fstream], [],
                 [AC_DEFINE([HAVE_NO_FSTREAM], 1, [<fstream> not available])],
                 [#include <cstdlib>])

AC_CHECK_HEADERS([sstream], [],
                 [AC_DEFINE([HAVE_NO_SSTREAM], 1, [<sstream> not available])],
                 [#include <cstdlib>])

TAC_ARG_WITH_AR

# ------------------------------------------------------------------------
# MPI link check
# ------------------------------------------------------------------------
TAC_ARG_CONFIG_MPI

TAC_ARG_CHECK_MPI

AM_CONDITIONAL(FEI_SERIAL, [test "X${HAVE_MPI}" != "Xyes"])

AC_CXX_STD_IOS_BASE_FMTFLAGS

AC_CXX_UNORDERED

# ------------------------------------------------------------------------
# Check for other Trilinos packages
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_CAN_USE_PACKAGE(fei, teuchos, TEUCHOS, BUILD_TEUCHOS, yes, [Enable Teuchos support.], [Enable Teuchos support. Can be overridden with --disable-fei-teuchos.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(fei, epetra, EPETRA, BUILD_EPETRA, yes, [Enable Epetra support.], [Enable Epetra support. Can be overridden with --disable-fei-epetra.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(fei, aztecoo, AZTECOO, BUILD_AZTECOO, yes, [Enable AztecOO support.], [Enable AztecOO support. Can be overridden with --disable-fei-aztecoo.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(fei, amesos, AMESOS, BUILD_AMESOS, yes, [Enable Amesos support.], [Enable Amesos support. Can be overridden with --disable-fei-amesos.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(fei, ml, ML, BUILD_ML, yes, [Enable ML support.], [Enable ML support. Can be overridden with --disable-fei-ml.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(fei, ifpack, IFPACK, BUILD_IFPACK, yes, [Enable Ifpack support.], [Enable Ifpack support. Can be overridden with --disable-fei-ifpack.])

# ------------------------------------------------------------------------
# Check for PETSc include directives
# ------------------------------------------------------------------------

SPECIFIED_PETSC_INCLUDE=no

AC_ARG_WITH(petsc-include,
  AC_HELP_STRING([--with-petsc-include=-IDIR], [Specify PETSc include]),
  [
  AC_MSG_CHECKING(whether petsc-include specified)
  AC_MSG_RESULT([${withval}])
  SPECIFIED_PETSC_INCLUDE=yes
  PETSCINC="${withval}"
  ]
)

AM_CONDITIONAL(FEI_HAVE_PETSC, [test "X${SPECIFIED_PETSC_INCLUDE}" != "Xno"])

if test "X${SPECIFIED_PETSC_INCLUDE}" != "Xno"; then
PETSC_INCLUDE=${PETSCINC}
AC_SUBST( PETSC_INCLUDE )
fi

# ------------------------------------------------------------------------
# Check for PETSC library directory
# ------------------------------------------------------------------------

SPECIFIED_PETSC_LIB=no

AC_ARG_WITH(petsc-libdir,
  AC_HELP_STRING([--with-petsc-libdir=LIBDIR], [Specify PETSc library directory]),
  [
  AC_MSG_CHECKING(whether petsc-libdir specified)
  AC_MSG_RESULT([${withval}])
  SPECIFIED_PETSC_LIB=yes
  PETSCLIBDIR="${withval}"
  ]
)

if test "X${SPECIFIED_PETSC_LIB}" != "Xno"; then
TAC_CHECK_FILE(${PETSCLIBDIR}/libpetsc.a)
PETSC_LIB_DIR=${PETSCLIBDIR}
AC_SUBST( PETSC_LIB_DIR )
AC_DEFINE([FEI_HAVE_PETSC],,
          [Enable use of PETSC library])
fi

# ------------------------------------------------------------------------
# Check for Prometheus include directives
# ------------------------------------------------------------------------

SPECIFIED_PROMETHEUS_INCLUDE=no

AC_ARG_WITH(prometheus-include,
  AC_HELP_STRING([--with-prometheus-include=-IDIR], [Specify Prometheus include]),
  [
  AC_MSG_CHECKING(whether prometheus-include specified)
  AC_MSG_RESULT([${withval}])
  SPECIFIED_PROMETHEUS_INCLUDE=yes
  PROMETHEUSINC="${withval}"
  ]
)

AM_CONDITIONAL(FEI_HAVE_PROMETHEUS, [test "X${SPECIFIED_PROMETHEUS_INCLUDE}" != "Xno"])

if test "X${SPECIFIED_PROMETHEUS_INCLUDE}" != "Xno"; then
PROMETHEUS_INCLUDE=${PROMETHEUSINC}
AC_SUBST( PROMETHEUS_INCLUDE )
fi

# ------------------------------------------------------------------------
# Check for Prometheus library directory
# ------------------------------------------------------------------------

SPECIFIED_PROMETHEUS_LIB=no

AC_ARG_WITH(prometheus-libdir,
  AC_HELP_STRING([--with-prometheus-libdir=LIBDIR], [Specify Prometheus library directory]),
  [
  AC_MSG_CHECKING(whether prometheus-libdir specified)
  AC_MSG_RESULT([${withval}])
  SPECIFIED_PROMETHEUS_LIB=yes
  PROMETHEUSLIBDIR="${withval}"
  ]
)

if test "X${SPECIFIED_PROMETHEUS_LIB}" != "Xno"; then
TAC_CHECK_FILE(${PROMETHEUSLIBDIR}/libprometheus.a)
PROMETHEUS_LIB_DIR=${PROMETHEUSLIBDIR}
AC_SUBST( PROMETHEUS_LIB_DIR )
AC_DEFINE([FEI_HAVE_PROMETHEUS],,
          [Enable use of PROMETHEUS library])
fi

# ------------------------------------------------------------------------
# Check for HYPRE include directives
# ------------------------------------------------------------------------

SPECIFIED_HYPRE_INCLUDE=no

AC_ARG_WITH(hypre-include,
  AC_HELP_STRING([--with-hypre-include=-IDIR], [Specify HYPRE include]),
  [
  AC_MSG_CHECKING(whether hypre-include specified)
  AC_MSG_RESULT([${withval}])
  SPECIFIED_HYPRE_INCLUDE=yes
  HYPREINC="${withval}"
  ]
)

AM_CONDITIONAL(FEI_HAVE_HYPRE, [test "X${SPECIFIED_HYPRE_INCLUDE}" != "Xno"])

if test "X${SPECIFIED_HYPRE_INCLUDE}" != "Xno"; then
HYPRE_INCLUDE=${HYPREINC}
AC_SUBST( HYPRE_INCLUDE )
fi

# ------------------------------------------------------------------------
# Check for HYPRE library directory
# ------------------------------------------------------------------------

SPECIFIED_HYPRE_LIB=no

AC_ARG_WITH(hypre-libdir,
  AC_HELP_STRING([--with-hypre-libdir=LIBDIR], [Specify HYPRE library directory]),
  [
  AC_MSG_CHECKING(whether hypre-libdir specified)
  AC_MSG_RESULT([${withval}])
  SPECIFIED_HYPRE_LIB=yes
  HYPRELIBDIR="${withval}"
  ]
)

if test "X${SPECIFIED_HYPRE_LIB}" != "Xno"; then
TAC_CHECK_FILE(${HYPRELIBDIR}/libHYPRE_FEI.a)
HYPRE_LIB_DIR=${HYPRELIBDIR}
AC_SUBST( HYPRE_LIB_DIR )
AC_DEFINE([FEI_HAVE_HYPRE],,
          [Enable use of HYPRE library])
fi

# ------------------------------------------------------------------------
# Check for FETI-DP include directives
# ------------------------------------------------------------------------

SPECIFIED_FETIDP_INCLUDE=no

AC_ARG_WITH(fetidp-include,
  AC_HELP_STRING([--with-fetidp-include=-IDIR], [Specify FETI-DP include]),
  [
  AC_MSG_CHECKING(whether fetidp-include specified)
  AC_MSG_RESULT([${withval}])
  SPECIFIED_FETIDP_INCLUDE=yes
  FETIDPINC="${withval}"
  ]
)

AM_CONDITIONAL(FEI_HAVE_FETI, [test "X${SPECIFIED_FETIDP_INCLUDE}" != "Xno"])

if test "X${SPECIFIED_FETIDP_INCLUDE}" != "Xno"; then
FETIDP_INCLUDE=${FETIDPINC}
AC_SUBST( FETIDP_INCLUDE )
fi

# ------------------------------------------------------------------------
# Check for FETI-DP library directory
# ------------------------------------------------------------------------

SPECIFIED_FETIDP_LIB=no

AC_ARG_WITH(fetidp-libdir,
  AC_HELP_STRING([--with-fetidp-libdir=LIBDIR], [Specify FETI-DP library directory]),
  [
  AC_MSG_CHECKING(whether fetidp-libdir specified)
  AC_MSG_RESULT([${withval}])
  SPECIFIED_FETIDP_LIB=yes
  FETIDPLIBDIR="${withval}"
  ]
)

if test "X${SPECIFIED_FETIDP_LIB}" != "Xno"; then
TAC_CHECK_FILE(${FETIDPLIBDIR}/libfeti-mpi.a)
FETIDP_LIB_DIR=${FETIDPLIBDIR}
AC_SUBST( FETIDP_LIB_DIR )
AC_DEFINE([FEI_HAVE_FETI],,
          [Enable use of FETI library])
fi

# ------------------------------------------------------------------------
# Check for --with-instrument-exe= command.
# ------------------------------------------------------------------------

SPECIFIED_INSTRUMENT=no

AC_ARG_WITH(instrument-exe,
  AC_HELP_STRING([--with-instrument-exe=CMD], [Specify cmd for instrumenting executable. Example: "purify -always-use-cache-dir -cache-dir=/path/dir"]),
  [
  AC_MSG_CHECKING(whether instrument-exe specified)
  AC_MSG_RESULT([${withval}])
  SPECIFIED_INSTRUMENT=yes
  INSTRUMENTCMD="${withval}"
  ]
)

AM_CONDITIONAL(FEI_INSTRUMENT, [test "X${SPECIFIED_INSTRUMENT}" != "Xno"])

if test "X${SPECIFIED_INSTRUMENT}" != "Xno"; then
INSTRUMENT_CMD=${INSTRUMENTCMD}
AC_SUBST( INSTRUMENT_CMD )
fi


# ------------------------------------------------------------------------
# Check for --with-run-cmd= command.
# ------------------------------------------------------------------------

SPEC_RUN_COMMAND=no

AC_ARG_WITH(run-cmd,
  AC_HELP_STRING([--with-run-cmd=CMD],
                 [Specify command for running executable (if not specified, defaults to "mpirun -np"). The form of the run-cmd should be such that a parallel program can be run using this syntax: "CMD <nprocs> <executable> <executable-args>"]),
  [
  AC_MSG_CHECKING(whether run-cmd specified)
  AC_MSG_RESULT([${withval}])
  SPEC_RUN_COMMAND=yes
  RUNCMD="${withval}"
  ]
)

AM_CONDITIONAL(RUN_COMMAND_SPECIFIED, [test "X${SPEC_RUN_COMMAND}" != "Xno"])

if test "X${SPEC_RUN_COMMAND}" != "Xno"; then
SPECIFIED_RUN_COMMAND=${RUNCMD}
AC_SUBST( SPECIFIED_RUN_COMMAND )
fi


# ------------------------------------------------------------------------
# Checks if tests and examples should be built
# ------------------------------------------------------------------------

#  These options can disable the tests and examples of a package.
TAC_ARG_ENABLE_FEATURE(tests, [Build tests for FEI package], TESTS, yes)
TAC_ARG_ENABLE_FEATURE(examples, [Build examples for FEI package], EXAMPLES, yes)

AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_tests" != "Xno")
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_examples" != "Xno")

#We now build tests and examples through separate make targets, rather than
#during "make".  We still need to conditionally include the test and example
#in SUBDIRS, even though SUB_TEST and SUB_EXAMPLE will never be
#defined, so that the tests and examples are included in the distribution
#tarball.
AM_CONDITIONAL(SUB_TEST, test "X$ac_cv_use_sub_test" = "Xyes")
AM_CONDITIONAL(SUB_EXAMPLE, test "X$ac_cv_use_sub_example" = "Xyes")

# ------------------------------------------------------------------------
# Specify other directories
# ------------------------------------------------------------------------

# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS


AC_SUBST(ac_aux_dir)

AC_CONFIG_FILES([
                Makefile
		Makefile.export.fei
                base/Makefile
                test_utils/Makefile
                support-Trilinos/Makefile
                support-PETSc/Makefile
                test/Makefile
                test/cr_lagrange_20quad/Makefile
                test/cr_lagrange_20quad_AZLSC/Makefile
                test/cr_lagrange_20quad_old/Makefile
                test/cr_multifield_slave/Makefile
                test/elemDOF_Aztec/Makefile
                test/remEssBCs/Makefile
		utest/Makefile
		utest_cases/Makefile
                example/Makefile
                example/beam/Makefile
                example/poisson/Makefile
                ])

AC_OUTPUT()

echo "----------------------------------------"
echo "Finished Running FEI Configure Script"
echo "----------------------------------------"

