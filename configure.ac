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

AM_INIT_AUTOMAKE(Trilinos, 2.1)

# Specify required version of autoconf.

AC_PREREQ(2.52)

# Must indicate for packages that can be configured from the package level that the 
# configure/build process is happening within Trilinos.

ac_configure_args="$ac_configure_args --with-trilinos"

TAC_ARG_ENABLE_FEATURE(referencey12m, [Build y12m lib included in the distribution], Y12M, no)
AM_CONDITIONAL(BUILD_Y12M, test "X$ac_cv_use_referencey12m" != "Xno")

# For now, blas and lapack will not be included
#TAC_ARG_ENABLE_FEATURE(referenceblasandlapack, [Build blas and lapack libs included in the distribution], BLASANDLAPACK, no)
#AM_CONDITIONAL(BUILD_BLAS, test "X$ac_cv_use_referenceblasandlapack" != "Xno")

#TAC_ARG_ENABLE_FEATURE(referencelapack, [Build lapack lib included in the distribution], LAPACK, no)
#AM_CONDITIONAL(BUILD_LAPACK, test "X$ac_cv_use_referencelapack" != "Xno" || test "X$ac_cv_use_referenceblasandlapack" != "Xno")

# After building Nox from the Trilinos level is stabilized, the default should be 
# changed to yes.
TAC_ARG_ENABLE_FEATURE(nox, [Build nox], NOX, no)
AM_CONDITIONAL(BUILD_NOX, test "X$ac_cv_use_nox" != "Xno")

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

TAC_ARG_ENABLE_FEATURE(epetra, [Build epetra], EPETRA, yes)
AM_CONDITIONAL(BUILD_EPETRA, test "X$ac_cv_use_epetra" != "Xno" || 
	test "X$ac_cv_use_ifpack" != "Xno" || 
	test "X$ac_cv_use_oldifpack" != "Xno" ||
	test "X$ac_cv_use_triutils" != "Xno" ||
	test "X$ac_cv_use_aztecoo" != "Xno" ||
	test "X$ac_cv_use_komplex" != "Xno")

if test "X$ac_cv_use_y12m" != "Xno"; then
  AC_CONFIG_SUBDIRS([external/y12m])
fi

# For now, blas and lapck will not be included
#if test "X$ac_cv_use_reference_blasandlapack" != "Xno"; then
#  AC_CONFIG_SUBDIRS([external/blas])
#fi

#if test "X$ac_cv_use_reference_lapack" != "Xno" || 
#	test "X$ac_cv_use_referenceblasandlapack" != "Xno"; then
#  AC_CONFIG_SUBDIRS([external/lapack])
#fi


if test "X$ac_cv_use_epetra" != "Xno" || test "X$ac_cv_use_ifpack" != "Xno" ||
	test "X$ac_cv_use_oldifpack" != "Xno" ||
	test "X$ac_cv_use_triutils" != "Xno" || 
	test "X$ac_cv_use_aztecoo" != "Xno" || 
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

AC_CONFIG_FILES([
		Makefile
		packages/Makefile
		external/Makefile])
 AC_OUTPUT