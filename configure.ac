#
#  Refer to packages/new_package/README for suggestions on 
#  adding a new package.  
#

# Process this file with autoconf to produce a configure script.
AC_INIT(Trilinos, 2.1, mheroux@cs.sandia.gov)
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

# Use automake, getting name and version from variables defined by
# AC_INIT command.

AM_INIT_AUTOMAKE(Trilinos, 3.0.1)

# Specify required version of autoconf.

AC_PREREQ(2.52)

# Must indicate for packages that can be configured from the package level that the 
# configure/build process is happening within Trilinos.

ac_configure_args="$ac_configure_args --with-trilinos"

TAC_ARG_ENABLE_FEATURE(referencey12m, [Build y12m lib included in the distribution], Y12M, no)
AM_CONDITIONAL(BUILD_Y12M, test "X$ac_cv_use_referencey12m" != "Xno")

TAC_ARG_ENABLE_FEATURE(amesos, [Build amesos], AMESOS, no)
AM_CONDITIONAL(BUILD_AMESOS, test "X$ac_cv_use_amesos" != "Xno")

TAC_ARG_ENABLE_FEATURE(nox, [Build nox], NOX, no)
AM_CONDITIONAL(BUILD_NOX, test "X$ac_cv_use_nox" != "Xno")

TAC_ARG_ENABLE_FEATURE(meros, [Build meros], MEROS, no)
AM_CONDITIONAL(BUILD_MEROS, test "X$ac_cv_use_meros" != "Xno")

TAC_ARG_ENABLE_FEATURE(tsf, [Build TSF], TSF, no)
AM_CONDITIONAL(BUILD_TSF, test "X$ac_cv_use_tsf" != "Xno" ||
	test "X$ac_cv_use_meros" != "Xno")

TAC_ARG_ENABLE_FEATURE(tsfcore, [Build TSFCore], TSFCORE, no)
AM_CONDITIONAL(BUILD_TSFCORE, test "X$ac_cv_use_tsfcore" != "Xno")

TAC_ARG_ENABLE_FEATURE(tsfcoreutils, [Build TSFCoreUtils], TSFCOREUTILS, no)
AM_CONDITIONAL(BUILD_TSFCOREUTILS, test "X$ac_cv_use_tsfcoreutils" != "Xno" ||
	test "X$ac_cv_use_tsfcore" != "Xno")

TAC_ARG_ENABLE_FEATURE(komplex, [Build komplex], KOMPLEX, yes)
AM_CONDITIONAL(BUILD_KOMPLEX, test "X$ac_cv_use_komplex" != "Xno")

TAC_ARG_ENABLE_FEATURE(aztecoo, [Build aztecoo], AZTECOO, yes)
AM_CONDITIONAL(BUILD_AZTECOO, test "X$ac_cv_use_aztecoo" != "Xno" || 
	test "X$ac_cv_use_komplex" != "Xno")

TAC_ARG_ENABLE_FEATURE(belos, [Build belos], BELOS, no)
AM_CONDITIONAL(BUILD_BELOS, test "X$ac_cv_use_belos" != "Xno")

TAC_ARG_ENABLE_FEATURE(anasazi, [Build anasazi], ANASAZI, no)
AM_CONDITIONAL(BUILD_ANASAZI, test "X$ac_cv_use_anasazi" != "Xno")

TAC_ARG_ENABLE_FEATURE(ifpack, [Build ifpack], IFPACK, yes)
TAC_ARG_ENABLE_FEATURE(oldifpack, [Build ifpack and old ifpack ], OLDIFPACK, no)
AM_CONDITIONAL(BUILD_IFPACK, test "X$ac_cv_use_ifpack" != "Xno" ||
	test "X$ac_cv_use_belos" != "Xno" ||
	test "X$ac_cv_use_oldifpack" != "Xno")

TAC_ARG_ENABLE_FEATURE(triutils, [Build triutils], TRIUTILS, yes)
AM_CONDITIONAL(BUILD_TRIUTILS, test "X$ac_cv_use_triutils" != "Xno" || 
	test "X$ac_cv_use_belos" != "Xno")

TAC_ARG_ENABLE_FEATURE(epetraext, [Build epetraext], EPETRAEXT, no)
AM_CONDITIONAL(BUILD_EPETRAEXT, test "X$ac_cv_use_epetraext" != "Xno")

## #np# - Copy the next two lines, and change the name
TAC_ARG_ENABLE_FEATURE(new_package, [Build new_package], NEW_PACKAGE, no)
AM_CONDITIONAL(BUILD_NEW_PACKAGE, test "X$ac_cv_use_new_package" != "Xno")



## #np# - If your package requires epetra - you should add
## it to the following list
##	test "X$ac_cv_use_new_package" != "Xno" || 
TAC_ARG_ENABLE_FEATURE(epetra, [Build epetra], EPETRA, yes)
AM_CONDITIONAL(BUILD_EPETRA, test "X$ac_cv_use_epetra" != "Xno" || 
	test "X$ac_cv_use_ifpack" != "Xno" || 
	test "X$ac_cv_use_tsf" != "Xno" || 
	test "X$ac_cv_use_amesos" != "Xno" || 
	test "X$ac_cv_use_oldifpack" != "Xno" ||
	test "X$ac_cv_use_triutils" != "Xno" ||
	test "X$ac_cv_use_aztecoo" != "Xno" ||
	test "X$ac_cv_use_komplex" != "Xno" ||
	test "X$ac_cv_use_belos" != "Xno" ||
	test "X$ac_cv_use_anasazi" != "Xno" ||
	test "X$ac_cv_use_new_package" != "Xno" || 
        test "X$ac_cv_use_epetraext" != "Xno"  ||
        test "X$ac_cv_use_meros" != "Xno")

if test "X$ac_cv_use_y12m" != "Xno"; then
  AC_CONFIG_SUBDIRS([external/y12m])
fi

if test "X$ac_cv_use_epetra" != "Xno" || test "X$ac_cv_use_ifpack" != "Xno" ||
	test "X$ac_cv_use_oldifpack" != "Xno" ||
	test "X$ac_cv_use_triutils" != "Xno" || 
	test "X$ac_cv_use_tsf" != "Xno" || 
	test "X$ac_cv_use_amesos" != "Xno" || 
	test "X$ac_cv_use_aztecoo" != "Xno" || 
	test "X$ac_cv_use_new_package" != "Xno" || 
	test "X$ac_cv_use_meros" != "Xno" || 
        test "X$ac_cv_use_epetraext" != "Xno" ||
	test "X$ac_cv_use_belos" != "Xno" ||
	test "X$ac_cv_use_komplex" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/epetra])
fi

if test "X$ac_cv_use_triutils" != "Xno" ||
	test "X$ac_cv_use_belos" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/triutils])
fi

if test "X$ac_cv_use_ifpack" != "Xno" ||
	test "X$ac_cv_use_belos" != "Xno" ||
	test "X$ac_cv_use_oldifpack" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/ifpack])
fi

if test "X$ac_cv_use_aztecoo" != "Xno" || 
	test "X$ac_cv_use_komplex" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/aztecoo])
fi

if test "X$ac_cv_use_komplex" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/komplex])
fi

if test "X$ac_cv_use_nox" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/nox])
fi

if test "X$ac_cv_use_tsf" != "Xno" ||
	test "X$ac_cv_use_meros" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/TSF])
fi

if test "X$ac_cv_use_tsfcore" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/TSFCore])
fi

if test "X$ac_cv_use_tsfcore" != "Xno" || 
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

if test "X$ac_cv_use_new_package" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/new_package])
fi

AC_ARG_ENABLE(,[

The options listed above are specific to the Trilinos level.
For a list of all options for all packages, type 
"configure --help=recursive".  The output for this command
is bulky.  Consider using the "more" or "grep" command to 
find a specific option.

To view a list of options for a specific package, cd to 
"Trilinos/package/package_name" (replacing "package_name"
with the name of the package) and type "configure --help".

A list of configuration options that are common to all
Trilinos packages can be found on the web at
"http://software.sandia.gov/Trilinos/config_options.html".
Many of the descriptions of the options found at this
webpage are more detailed than those found using
"configure --help".])

AC_CONFIG_FILES([
		Makefile
		packages/Makefile
		external/Makefile])
 AC_OUTPUT
