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

# After building Nox from the Trilinos level is stabilized, the default should be 
# changed to yes.
TAC_ARG_ENABLE_FEATURE(nox, [Build nox], NOX, no)
AM_CONDITIONAL(BUILD_NOX, test "X$ac_cv_use_nox" != "Xno")

# After building TSF from the Trilinos level is stabilized, 
# the default should be changed to yes.
TAC_ARG_ENABLE_FEATURE(tsf, [Build TSF], TSF, no)
AM_CONDITIONAL(BUILD_TSF, test "X$ac_cv_use_tsf" != "Xno")

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

TAC_ARG_ENABLE_FEATURE(ifpack, [Build ifpack], IFPACK, yes)
TAC_ARG_ENABLE_FEATURE(oldifpack, [Build ifpack and old ifpack ], OLDIFPACK, no)
AM_CONDITIONAL(BUILD_IFPACK, test "X$ac_cv_use_ifpack" != "Xno" ||
	test "X$ac_cv_use_oldifpack" != "Xno")

TAC_ARG_ENABLE_FEATURE(triutils, [Build triutils], TRIUTILS, yes)
AM_CONDITIONAL(BUILD_TRIUTILS, test "X$ac_cv_use_triutils" != "Xno")

# After building EpetraExt from the Trilinos level is stabilized, the 
#default could be changed to yes.
TAC_ARG_ENABLE_FEATURE(epetraext, [Build epetraext], EPETRAEXT, no)
AM_CONDITIONAL(BUILD_EPETRAEXT, test "X$ac_cv_use_epetraext" != "Xno")

## #np# - Copy the next four lines, and change the name
## # After building New_Package from the Trilinos level is stabilized, the 
## #default could be changed to yes.
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
	test "X$ac_cv_use_new_package" != "Xno" || 
        test "X$ac_cv_use_epetraext" != "Xno")

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
        test "X$ac_cv_use_epetraext" != "Xno" ||
	test "X$ac_cv_use_komplex" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/epetra])
fi

if test "X$ac_cv_use_triutils" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/triutils])
fi

if test "X$ac_cv_use_ifpack" != "Xno" ||
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

if test "X$ac_cv_use_tsf" != "Xno"; then
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

# REQUIRED CHANGE 
# Copy the next three lines change the name
if test "X$ac_cv_use_new_package" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/new_package])
fi

AC_CONFIG_FILES([
		Makefile
		packages/Makefile
		external/Makefile])
 AC_OUTPUT
