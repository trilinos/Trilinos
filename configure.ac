# @HEADER
# @HEADER

--#
#  Refer to packages/new_package/README for suggestions on 
#  adding a new package.  
#

# Process this file with autoconf to produce a configure script.
AC_INIT(Trilinos, 5.0d, maherou@sandia.gov)
AC_CONFIG_SRCDIR([packages/epetra/src/Epetra_Object.h])
AM_CONFIG_HEADER([Trilinos_config.h:Trilinos_config.h.in])

# AM_MAINTAINER_MODE turns off maintainer-only makefile targets by
# default, and changes configure to understand a
# --enable-maintainer-mode option. --enable-maintainer-mode turns the
# maintainer-only targets back on. The maintainer-only makefile
# targets permit end users to clean automatically-generated files such
# as configure, which means they have to have autoconf and automake
# installed to repair the damage. AM_MAINTAINER_MODE makes it a bit
# harder for users to shoot themselves in the foot.

AM_MAINTAINER_MODE

# Specify required version of automake.

AM_INIT_AUTOMAKE(1.8.2)

# Specify required version of autoconf.

AC_PREREQ(2.59)

# Must indicate for packages that can be configured from the package level that the 
# configure/build process is happening within Trilinos.

ac_configure_args="$ac_configure_args --with-trilinos"

echo "NOTE: If a package is not enabled based on the"
echo "configure options, but is required by another"
echo "package that has been enabled, it will"
echo "automatically be configured and built, even if"
echo "the list below does not indicate so."

TAC_ARG_ENABLE_FEATURE_SUB( default, packages, [Build all default packages unless explicitly disabled.  When this option is disabled, only packages explicitly enabled and packages needed to satisfy the dependencies of enabled packages will be built], DEFAULT_PACKAGES, yes)
 
TAC_ARG_ENABLE_FEATURE(amesos, [Build amesos], AMESOS, no)
AM_CONDITIONAL(BUILD_AMESOS, test "X$ac_cv_use_amesos" != "Xno")

TAC_ARG_ENABLE_DEFAULT_FEATURE(loca, [Build loca], LOCA)
#The following line is not necessary for LOCA because LOCA
#lives in the NOX directory structure
#AM_CONDITIONAL(BUILD_NOX, test "X$ac_cv_use_nox" != "Xno")

TAC_ARG_ENABLE_DEFAULT_FEATURE(nox, [Build nox], NOX)
AM_CONDITIONAL(BUILD_NOX, test "X$ac_cv_use_nox" != "Xno" ||
	test "X$ac_cv_use_loca" != "Xno")

TAC_ARG_ENABLE_FEATURE(meros, [Build meros], MEROS, no)
AM_CONDITIONAL(BUILD_MEROS, test "X$ac_cv_use_meros" != "Xno")

TAC_ARG_ENABLE_FEATURE(tsf, [Build TSF], TSF, no)
AM_CONDITIONAL(BUILD_TSF, test "X$ac_cv_use_tsf" != "Xno" ||
	test "X$ac_cv_use_meros" != "Xno")

TAC_ARG_ENABLE_FEATURE(tsfextended, [Build TSFExtended], TSFEXTENDED, no)
AM_CONDITIONAL(BUILD_TSFEXTENDED, test "X$ac_cv_use_tsfextended" != "Xno")

TAC_ARG_ENABLE_FEATURE(tsfcore, [Build TSFCore], TSFCORE, no)
AM_CONDITIONAL(BUILD_TSFCORE, test "X$ac_cv_use_tsfcore" != "Xno" || 
  test "X$ac_cv_use_tsfextended" != "Xno")

TAC_ARG_ENABLE_FEATURE(tsfcoreutils, [Build TSFCoreUtils], TSFCOREUTILS, no)
AM_CONDITIONAL(BUILD_TSFCOREUTILS, test "X$ac_cv_use_tsfcoreutils" != "Xno" ||
	test "X$ac_cv_use_tsfcore" != "Xno" || 
  test "X$ac_cv_use_tsfextended" != "Xno")

TAC_ARG_ENABLE_FEATURE(komplex, [Build komplex], KOMPLEX, no)
AM_CONDITIONAL(BUILD_KOMPLEX, test "X$ac_cv_use_komplex" != "Xno")

TAC_ARG_ENABLE_DEFAULT_FEATURE(aztecoo, [Build aztecoo], AZTECOO)
AM_CONDITIONAL(BUILD_AZTECOO, test "X$ac_cv_use_aztecoo" != "Xno" || 
	test "X$ac_cv_use_komplex" != "Xno")

TAC_ARG_ENABLE_FEATURE(belos, [Build belos], BELOS, no)
AM_CONDITIONAL(BUILD_BELOS, test "X$ac_cv_use_belos" != "Xno")

TAC_ARG_ENABLE_FEATURE(anasazi, [Build anasazi], ANASAZI, no)
AM_CONDITIONAL(BUILD_ANASAZI, test "X$ac_cv_use_anasazi" != "Xno")

TAC_ARG_ENABLE_DEFAULT_FEATURE(ifpack, [Build ifpack], IFPACK)
AM_CONDITIONAL(BUILD_IFPACK, test "X$ac_cv_use_ifpack" != "Xno")

TAC_ARG_ENABLE_FEATURE(triutils, [Build triutils], TRIUTILS, no)
AM_CONDITIONAL(BUILD_TRIUTILS, test "X$ac_cv_use_triutils" != "Xno" || 
	test "X$ac_cv_use_aztecoo" != "Xno" || 
	test "X$ac_cv_use_anasazi" != "Xno" || 
	test "X$ac_cv_use_belos" != "Xno" || 
	test "X$ac_cv_use_amesos" != "Xno")

TAC_ARG_ENABLE_DEFAULT_FEATURE(epetraext, [Build epetraext], EPETRAEXT)
AM_CONDITIONAL(BUILD_EPETRAEXT, test "X$ac_cv_use_epetraext" != "Xno")

TAC_ARG_ENABLE_FEATURE(tpetra, [Build tpetra], TPETRA, no)
AM_CONDITIONAL(BUILD_TPETRA, test "X$ac_cv_use_tpetra" != "Xno")

TAC_ARG_ENABLE_FEATURE(teuchos, [Build teuchos], TEUCHOS, no)
AM_CONDITIONAL(BUILD_TEUCHOS, test "X$ac_cv_use_teuchos" != "Xno" ||
	test "X$ac_cv_use_tpetra" != "Xno" || 
	test "X$ac_cv_use_tsfcore" != "Xno" || 
	test "X$ac_cv_use_tsfcoreutils" != "Xno" || 
	test "X$ac_cv_use_tsfextended" != "Xno" || 
        test "X$ac_cv_use_loca" != "Xno" ||
 	test "X$ac_cv_use_belos" != "Xno" || 
	test "X$ac_cv_use_amesos" != "Xno" || 
	test "X$ac_cv_use_anasazi" != "Xno") 

TAC_ARG_ENABLE_DEFAULT_FEATURE(ml, [Build ml], ML)
AM_CONDITIONAL(BUILD_ML, test "X$ac_cv_use_ml" != "Xno")

## #np# - Copy the next two lines, and change the name
TAC_ARG_ENABLE_DEFAULT_FEATURE(new_package, [Build new_package], NEW_PACKAGE)
AM_CONDITIONAL(BUILD_NEW_PACKAGE, test "X$ac_cv_use_new_package" != "Xno")


if test "X$ac_cv_use_teuchos" != "Xno" ||
	test "X$ac_cv_use_tpetra" != "Xno" ||
        test "X$ac_cv_use_tsfextended" != "Xno" ||
 	test "X$ac_cv_use_tsfcore" != "Xno" ||
	test "X$ac_cv_use_tsfcoreutils" != "Xno" ||
        test "X$ac_cv_use_loca" != "Xno" ||
	test "X$ac_cv_use_amesos" != "Xno" || 
 	test "X$ac_cv_use_anasazi" != "Xno" ||
	test "X$ac_cv_use_belos" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/teuchos])
fi

## #np# - If your package requires epetra - you should add
## it to the following list
##	test "X$ac_cv_use_new_package" != "Xno" || 
TAC_ARG_ENABLE_DEFAULT_FEATURE(epetra, [Build epetra], EPETRA)
AM_CONDITIONAL(BUILD_EPETRA, test "X$ac_cv_use_epetra" != "Xno" || 
	test "X$ac_cv_use_ifpack" != "Xno" || 
	test "X$ac_cv_use_tsf" != "Xno" || 
	test "X$ac_cv_use_tsfextended" != "Xno" || 
	test "X$ac_cv_use_amesos" != "Xno" || 
	test "X$ac_cv_use_triutils" != "Xno" ||
	test "X$ac_cv_use_aztecoo" != "Xno" ||
	test "X$ac_cv_use_komplex" != "Xno" ||
	test "X$ac_cv_use_belos" != "Xno" ||
	test "X$ac_cv_use_anasazi" != "Xno" ||
	test "X$ac_cv_use_new_package" != "Xno" || 
        test "X$ac_cv_use_epetraext" != "Xno"  ||
        test "X$ac_cv_use_meros" != "Xno")

TAC_ARG_ENABLE_FEATURE(kokkos, [Build kokkos], KOKKOS, no)
AM_CONDITIONAL(BUILD_KOKKOS, test "X$ac_cv_use_kokkos" != "Xno" ||
	test "X$ac_cv_use_tpetra" != "Xno")

if test "X$ac_cv_use_kokkos" != "Xno" ||
  	test "X$ac_cv_use_tpetra" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/kokkos])
fi

if test "X$ac_cv_use_epetra" != "Xno" || test "X$ac_cv_use_ifpack" != "Xno" ||
	test "X$ac_cv_use_triutils" != "Xno" || 
	test "X$ac_cv_use_tsf" != "Xno" || 
        test "X$ac_cv_use_tsfextended" != "Xno" ||
 	test "X$ac_cv_use_amesos" != "Xno" || 
	test "X$ac_cv_use_aztecoo" != "Xno" || 
	test "X$ac_cv_use_new_package" != "Xno" || 
	test "X$ac_cv_use_meros" != "Xno" || 
        test "X$ac_cv_use_epetraext" != "Xno" ||
	test "X$ac_cv_use_belos" != "Xno" ||
	test "X$ac_cv_use_anasazi" != "Xno" ||
	test "X$ac_cv_use_komplex" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/epetra])
fi

if test "X$ac_cv_use_triutils" != "Xno" ||
	test "X$ac_cv_use_aztecoo" != "Xno" ||
	test "X$ac_cv_use_anasazi" != "Xno" ||
	test "X$ac_cv_use_belos" != "Xno" ||
	test "X$ac_cv_use_amesos" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/triutils])
fi

if test "X$ac_cv_use_ifpack" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/ifpack])
fi

if test "X$ac_cv_use_aztecoo" != "Xno" || 
	test "X$ac_cv_use_komplex" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/aztecoo])
fi

if test "X$ac_cv_use_komplex" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/komplex])
fi

if test "X$ac_cv_use_tsf" != "Xno" ||
	test "X$ac_cv_use_meros" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/TSF])
fi

if test "X$ac_cv_use_tsfcore" != "Xno" ||
	test "X$ac_cv_use_tsfextended" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/TSFCore])
fi

if test "X$ac_cv_use_tsfextended" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/TSFExtended])
fi

if test "X$ac_cv_use_tsfcore" != "Xno" || 
        test "X$ac_cv_use_tsfextended" != "Xno" ||
	test "X$ac_cv_use_tsfcoreutils" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/TSFCoreUtils])
fi

if test "X$ac_cv_use_amesos" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/amesos])
fi

if test "X$ac_cv_use_epetraext" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/epetraext])
fi

if test "X$ac_cv_use_belos" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/belos])
fi

if test "X$ac_cv_use_anasazi" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/anasazi])
fi

if test "X$ac_cv_use_meros" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/meros])
fi

if test "X$ac_cv_use_ml" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/ml])
fi

if test "X$ac_cv_use_nox" != "Xno" ||
        test "X$ac_cv_use_loca" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/nox])
fi

## #np# - Copy the next three lines, and change package name
if test "X$ac_cv_use_new_package" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/new_package])
fi

if test "X$ac_cv_use_tpetra" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/tpetra])
fi

AC_ARG_ENABLE(,[

CAN'T FIND THE OPTION YOU NEED?

The options listed above are specific to the Trilinos level.
These options are limited to deciding which packages should
be built.  All other options such as compilers and compiler
flags are set at the package level.
For a list of all options for all packages, type 
"configure --help=recursive".  The output for this command
is bulky.  Consider using the "more" or "grep" command to 
find a specific option.

To view a list of options for a specific package, cd to 
"Trilinos/package/package_name" (replacing "package_name"
with the name of the package) and type "configure --help".

A list of configuration options that are common to all
Trilinos packages can be found on the web at
"http://software.sandia.gov/trilinos/configuration_manual.html".
Many of the descriptions of the options found at this
webpage are more detailed than those found using
"configure --help".

THE TRILINOS USER GUIDE contains a lot of additional information
about configuring and build Trilinos, as well as other
topics that will be helpful for Trilinos Users.  The Trilinos
User Guide can be found at:
"http://software.sandia.gov/trilinos/publications.html". 
])
AC_CONFIG_FILES([
		Makefile
		packages/Makefile])
 AC_OUTPUT
