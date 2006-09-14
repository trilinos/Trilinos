# @HEADER
# @HEADER

#
#  Refer to packages/new_package/README for suggestions on 
#  adding a new package.  
#

# Process this file with autoconf to produce a configure script.
AC_INIT(Trilinos, 8.0d, maherou@sandia.gov)
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

TAC_ARG_ENABLE_FEATURE(anasazi, [Build anasazi], ANASAZI, no)

TAC_ARG_ENABLE_DEFAULT_FEATURE(aztecoo, [Build aztecoo], AZTECOO)

TAC_ARG_ENABLE_FEATURE(belos, [Build belos], BELOS, no)

TAC_ARG_ENABLE_FEATURE(capo, [Build capo], CAPO, no)

TAC_ARG_ENABLE_FEATURE(claps, [Build claps], CLAPS, no)

TAC_ARG_ENABLE_FEATURE(didasko, [Build didasko (a Trilinos tutorial)], DIDASKO, no)

TAC_ARG_ENABLE_DEFAULT_FEATURE(epetra, [Build epetra], EPETRA)

TAC_ARG_ENABLE_DEFAULT_FEATURE(epetraext, [Build epetraext], EPETRAEXT)

#The external package is a mechanism for adding external packages to the Trilinos
#build process.
TAC_ARG_ENABLE_FEATURE(external, [Build external package in packages/external], EXTERNAL, no)

TAC_ARG_ENABLE_DEFAULT_FEATURE(ifpack, [Build ifpack], IFPACK)

TAC_ARG_ENABLE_FEATURE(isorropia, [Build isorropia], ISORROPIA, no)

TAC_ARG_ENABLE_FEATURE(kokkos, [Build kokkos], KOKKOS, no)

TAC_ARG_ENABLE_FEATURE(komplex, [Build komplex], KOMPLEX, no)

TAC_ARG_ENABLE_DEFAULT_FEATURE(loca, [Build loca], LOCA)

TAC_ARG_ENABLE_FEATURE(meros, [Build meros], MEROS, no)

TAC_ARG_ENABLE_FEATURE(moertel, [Build moertel], MOERTEL, no)

TAC_ARG_ENABLE_DEFAULT_FEATURE(ml, [Build ml], ML)

## #np# - Copy the next line, and change the name.  Place the line in 
## #np# - alphabetical order based on package name.
TAC_ARG_ENABLE_DEFAULT_FEATURE(new_package, [Build new_package], NEW_PACKAGE)

TAC_ARG_ENABLE_DEFAULT_FEATURE(nox, [Build nox], NOX)

TAC_ARG_ENABLE_FEATURE(pliris, [Build pliris], PLIRIS, no)

TAC_ARG_ENABLE_FEATURE(pytrilinos, [Build PyTrilinos], PYTRILINOS, no)

TAC_ARG_ENABLE_FEATURE(python, [Build PyTrilinos], EPYTHON, no)

if test "X$ac_cv_use_python" != "Xno"; then
  ac_hold_python="yes"
fi

TAC_ARG_WITH_PACKAGE(python, [Build PyTrilinos], WPYTHON, no)

#We need to do this to support with and enable options with the same argument.
#Eventually, we should maybe change the base m4's if this situation becomes
#common.

if test "X$ac_cv_use_python" != "Xyes" && test "X$ac_hold_python" != "Xyes"; then
  ac_cv_use_python="no"
else
  ac_cv_use_python="yes"
fi

TAC_ARG_ENABLE_FEATURE(rythmos, [Build rythmos], RYTHMOS, no)

TAC_ARG_ENABLE_FEATURE(sacado, [Build sacado], SACADO, no)

TAC_ARG_ENABLE_FEATURE(teuchos, [Build teuchos], TEUCHOS, no)

TAC_ARG_ENABLE_FEATURE(rtop, [Build RTOP], RTOP, no)

TAC_ARG_ENABLE_FEATURE(thyra, [Build thyra], THYRA, no)

TAC_ARG_ENABLE_FEATURE(stratimikos, [Build stratimikos], STRATIMIKOS, no)

TAC_ARG_ENABLE_FEATURE(moocho, [Build moocho], MOOCHO, no)

TAC_ARG_ENABLE_FEATURE(tpetra, [Build tpetra], TPETRA, no)

TAC_ARG_ENABLE_FEATURE(triutils, [Build triutils], TRIUTILS, no)

TAC_ARG_ENABLE_FEATURE(galeri, [Build galeri], GALERI, no)

TAC_ARG_ENABLE_FEATURE(webtrilinos, [Build WebTrilinos], WEBTRILINOS, no)

AC_ARG_ENABLE(aztecoo-teuchos,
 [AC_HELP_STRING([--enable-aztecoo-teuchos],[Enable AztecOO Teuchos::ParameterList interface. Default is no. Requires Teuchos package.])],
 ac_cv_use_aztecoo_teuchos=$enableval, ac_cv_use_aztecoo_teuchos=no)

if test "X$ac_cv_use_teuchos" != "Xno" ||
	test "X$ac_cv_use_thyra" != "Xno" ||
	test "X$ac_cv_use_stratimikos" != "Xno" ||
	test "X$ac_cv_use_tpetra" != "Xno" ||
        test "X$ac_cv_use_loca" != "Xno" ||
	test "X$ac_cv_use_galeri" != "Xno" || 
	test "X$ac_cv_use_amesos" != "Xno" || 
 	test "X$ac_cv_use_anasazi" != "Xno" ||
	test "X$ac_cv_use_aztecoo_teuchos" != "Xno" ||
	test "X$ac_cv_use_capo" != "Xno" ||
	test "X$ac_cv_use_epetraext" != "Xno" ||
	test "X$ac_cv_use_meros" != "Xno" ||
	test "X$ac_cv_use_moocho" != "Xno" ||
	test "X$ac_cv_use_nox" != "Xno" ||
	test "X$ac_cv_use_rtop" != "Xno" ||
	test "X$ac_cv_use_rythmos" != "Xno" ||
	test "X$ac_cv_use_isorropia" != "Xno" ||
	test "X$ac_cv_use_moertel" != "Xno" ||
	test "X$ac_cv_use_webtrilinos" != "Xno" ||
	test "X$ac_cv_use_belos" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/teuchos])
  BUILD_TEUCHOS=yes
  if test "X$ac_cv_use_teuchos_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-teuchos"
  fi
else
  BUILD_TEUCHOS=no
fi
AM_CONDITIONAL(BUILD_TEUCHOS, test "X$BUILD_TEUCHOS" != "Xno")

if test "X$ac_cv_use_rtop" != "Xno" ||
	test "X$ac_cv_use_moocho" != "Xno" ||
	test "X$ac_cv_use_capo" != "Xno" ||
	test "X$ac_cv_use_meros" != "Xno" ||
	test "X$ac_cv_use_rythmos" != "Xno" ||
	test "X$ac_cv_use_stratimikos" != "Xno" ||
	test "X$ac_cv_use_thyra" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/rtop])
  BUILD_RTOP=yes
  if test "X$ac_cv_use_rtop_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-rtop"
  fi
else
  BUILD_RTOP=no
fi
AM_CONDITIONAL(BUILD_RTOP, test "X$BUILD_RTOP" != "Xno")

if test "X$ac_cv_use_thyra" != "Xno" ||
	test "X$ac_cv_use_capo" != "Xno" ||
	test "X$ac_cv_use_meros" != "Xno" ||
	test "X$ac_cv_use_rythmos" != "Xno" ||
	test "X$ac_cv_use_stratimikos" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/thyra])
  BUILD_THYRA=yes
  if test "X$ac_cv_use_thyra_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-thyra"
  fi
else
  BUILD_THYRA=no
fi
AM_CONDITIONAL(BUILD_THYRA, test "X$BUILD_THYRA" != "Xno")

if test "X$ac_cv_use_stratimikos" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/stratimikos])
  BUILD_STRATIMIKOS=yes
  if test "X$ac_cv_use_stratimikos_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-stratimikos"
  fi
else
  BUILD_STRATIMIKOS=no
fi
AM_CONDITIONAL(BUILD_STRATIMIKOS, test "X$BUILD_STRATIMIKOS" != "Xno")

if test "X$ac_cv_use_kokkos" != "Xno" ||
  	test "X$ac_cv_use_tpetra" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/kokkos])
  BUILD_KOKKOS=yes
  if test "X$ac_cv_use_kokkos_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-kokkos"
  fi
else
  BUILD_KOKKOS=no
fi
AM_CONDITIONAL(BUILD_KOKKOS, test "X$BUILD_KOKKOS" != "Xno")

if test "X$ac_cv_use_epetra" != "Xno" || test "X$ac_cv_use_ifpack" != "Xno" ||
	test "X$ac_cv_use_triutils" != "Xno" || 
	test "X$ac_cv_use_galeri" != "Xno" || 
 	test "X$ac_cv_use_amesos" != "Xno" || 
	test "X$ac_cv_use_aztecoo" != "Xno" || 
	test "X$ac_cv_use_new_package" != "Xno" || 
	test "X$ac_cv_use_meros" != "Xno" || 
        test "X$ac_cv_use_epetraext" != "Xno" ||
	test "X$ac_cv_use_belos" != "Xno" ||
	test "X$ac_cv_use_anasazi" != "Xno" ||
        test "X$ac_cv_use_claps" != "Xno" ||
	test "X$ac_cv_use_meros" != "Xno" ||
        test "X$ac_cv_use_moertel" != "Xno" ||
	test "X$ac_cv_use_webtrilinos" != "Xno" ||
	test "X$ac_cv_use_komplex" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/epetra])
  BUILD_EPETRA=yes
  if test "X$ac_cv_use_epetra_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-epetra"
  fi
else
  BUILD_EPETRA=no
fi
AM_CONDITIONAL(BUILD_EPETRA, test "X$BUILD_EPETRA" != "Xno")

if test "X$ac_cv_use_pliris" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/pliris])
  BUILD_PLIRIS=yes
  if test "X$ac_cv_use_pliris_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-pliris"
  fi
else
  BUILD_PLIRIS=no
fi
AM_CONDITIONAL(BUILD_PLIRIS, test "X$BUILD_PLIRIS" != "Xno")
                                                                                                    
if test "X$ac_cv_use_triutils" != "Xno" ||
	test "X$ac_cv_use_moertel" != "Xno" ||
	test "X$ac_cv_use_komplex" != "Xno" ||
	test "X$ac_cv_use_webtrilinos" != "Xno" ||
	test "X$ac_cv_use_meros" != "Xno" ||
	test "X$ac_cv_use_aztecoo" != "Xno"; then 
  AC_CONFIG_SUBDIRS([packages/triutils])
  BUILD_TRIUTILS=yes
  if test "X$ac_cv_use_triutils_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-triutils"
  fi
else
  BUILD_TRIUTILS=no
fi
AM_CONDITIONAL(BUILD_TRIUTILS, test "X$BUILD_TRIUTILS" != "Xno")

if test "X$ac_cv_use_galeri" != "Xno" ||
	test "X$ac_cv_use_webtrilinos" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/galeri])
  BUILD_GALERI=yes
  if test "X$ac_cv_use_galeri_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-galeri"
  fi
else
  BUILD_GALERI=no
fi
AM_CONDITIONAL(BUILD_GALERI, test "X$BUILD_GALERI" != "Xno")

if test "X$ac_cv_use_ifpack" != "Xno" ||
	test "X$ac_cv_use_webtrilinos" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/ifpack])
  BUILD_IFPACK=yes
  if test "X$ac_cv_use_ifpack_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-ifpack"
  fi
else
  BUILD_IFPACK=no
fi
AM_CONDITIONAL(BUILD_IFPACK, test "X$BUILD_IFPACK" != "Xno")

if test "X$ac_cv_use_aztecoo" != "Xno" || 
	test "X$ac_cv_use_moertel" != "Xno" ||
	test "X$ac_cv_use_meros" != "Xno" ||
	test "X$ac_cv_use_webtrilinos" != "Xno" ||
	test "X$ac_cv_use_komplex" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/aztecoo])
  BUILD_AZTECOO=yes
  if test "X$ac_cv_use_aztecoo_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-aztecoo"
  fi
else
  BUILD_AZTECOO=no
fi
AM_CONDITIONAL(BUILD_AZTECOO, test "X$BUILD_AZTECOO" != "Xno")

if test "X$ac_cv_use_claps" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/claps])
  BUILD_CLAPS=yes
  if test "X$ac_cv_use_claps_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-claps"
  fi
else
  BUILD_CLAPS=no
fi
AM_CONDITIONAL(BUILD_CLAPS, test "X$BUILD_CLAPS" != "Xno")

if test "X$ac_cv_use_komplex" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/komplex])
  BUILD_KOMPLEX=yes
  if test "X$ac_cv_use_komplex_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-komplex"
  fi
else
  BUILD_KOMPLEX=no
fi
AM_CONDITIONAL(BUILD_KOMPLEX, test "X$BUILD_KOMPLEX" != "Xno")

if test "X$ac_cv_use_amesos" != "Xno" ||
	test "X$ac_cv_use_webtrilinos" != "Xno" ||
        test "X$ac_cv_use_moertel" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/amesos])
  BUILD_AMESOS=yes
  if test "X$ac_cv_use_amesos_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-amesos"
  fi
else
  BUILD_AMESOS=no
fi
AM_CONDITIONAL(BUILD_AMESOS, test "X$BUILD_AMESOS" != "Xno")

if test "X$ac_cv_use_epetraext" != "Xno" ||
        test "X$ac_cv_use_moertel" != "Xno" ||
	test "X$ac_cv_use_webtrilinos" != "Xno" ||
	( test "X$ac_cv_use_thyra" != "Xno" && test "X$ac_cv_use_aztecoo" != "Xno" ); then
  AC_CONFIG_SUBDIRS([packages/epetraext])
  BUILD_EPETRAEXT=yes
  if test "X$ac_cv_use_epetraext_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-epetraext"
  fi
else
  BUILD_EPETRAEXT=no
fi
AM_CONDITIONAL(BUILD_EPETRAEXT, test "X$BUILD_EPETRAEXT" != "Xno")

if test "X$ac_cv_use_belos" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/belos])
  BUILD_BELOS=yes
  if test "X$ac_cv_use_belos_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-belos"
  fi
else
  BUILD_BELOS=no
fi
AM_CONDITIONAL(BUILD_BELOS, test "X$BUILD_BELOS" != "Xno")

if test "X$ac_cv_use_anasazi" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/anasazi])
  BUILD_ANASAZI=yes
  if test "X$ac_cv_use_anasazi_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-anasazi"
  fi
else
  BUILD_ANASAZI=no
fi
AM_CONDITIONAL(BUILD_ANASAZI, test "X$BUILD_ANASAZI" != "Xno")

if test "X$ac_cv_use_meros" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/meros])
  BUILD_MEROS=yes
  if test "X$ac_cv_use_meros_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-meros"
  fi
else
  BUILD_MEROS=no
fi
AM_CONDITIONAL(BUILD_MEROS, test "X$BUILD_MEROS" != "Xno")

if test "X$ac_cv_use_ml" != "Xno" ||
	test "X$ac_cv_use_webtrilinos" != "Xno" ||
	test "X$ac_cv_use_meros" != "Xno" ||
        test "X$ac_cv_use_moertel" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/ml])
  BUILD_ML=yes
  if test "X$ac_cv_use_ml_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-ml"
  fi
else
  BUILD_ML=no
fi
AM_CONDITIONAL(BUILD_ML, test "X$BUILD_ML" != "Xno")

if test "X$ac_cv_use_nox" != "Xno" ||
        test "X$ac_cv_use_loca" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/nox])
  BUILD_NOX=yes
  if test "X$ac_cv_use_nox_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-nox"
  fi
else
  BUILD_NOX=no
fi
AM_CONDITIONAL(BUILD_NOX, test "X$BUILD_NOX" != "Xno")

if test "X$ac_cv_use_isorropia" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/isorropia])
  BUILD_ISORROPIA=yes
  if test "X$ac_cv_use_isorropia_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-isorropia"
  fi
else
  BUILD_ISORROPIA=no
fi
AM_CONDITIONAL(BUILD_ISORROPIA, test "X$BUILD_ISORROPIA" != "Xno")

## #np# - Copy the next section, and change package name
if test "X$ac_cv_use_new_package" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/new_package])
  BUILD_NEW_PACKAGE=yes
  if test "X$ac_cv_use_new_package_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-new_package"
  fi
else
  BUILD_NEW_PACKAGE=no
fi
AM_CONDITIONAL(BUILD_NEW_PACKAGE, test "X$BUILD_NEW_PACKAGE" != "Xno")

if test "X$ac_cv_use_capo" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/capo])
  BUILD_CAPO=yes
  if test "X$ac_cv_use_capo_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-capo"
  fi
else
  BUILD_CAPO=no
fi
AM_CONDITIONAL(BUILD_CAPO, test "X$BUILD_CAPO" != "Xno")

if test "X$ac_cv_use_moertel" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/moertel])
  BUILD_MOERTEL=yes
  if test "X$ac_cv_use_moertel_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-moertel"
  fi
else
  BUILD_MOERTEL=no
fi
AM_CONDITIONAL(BUILD_MOERTEL, test "X$BUILD_MOERTEL" != "Xno")

if test "X$ac_cv_use_pytrilinos" != "Xno" ||
	test "X$ac_cv_use_webtrilinos" != "Xno" ||
	test "X$ac_cv_use_python" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/PyTrilinos])
  BUILD_PYTRILINOS=yes
  if test "X$ac_cv_use_pytrilinos_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-pytrilinos"
  fi
else
  BUILD_PYTRILINOS=no
fi
AM_CONDITIONAL(BUILD_PYTRILINOS, test "X$BUILD_PYTRILINOS" != "Xno")

if test "X$ac_cv_use_rythmos" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/rythmos])
  BUILD_RYTHMOS=yes
  if test "X$ac_cv_use_rythmos_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-rythmos"
  fi
else
  BUILD_RYTHMOS=no
fi
AM_CONDITIONAL(BUILD_RYTHMOS, test "X$BUILD_RYTHMOS" != "Xno")

if test "X$ac_cv_use_tpetra" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/tpetra])
  BUILD_TPETRA=yes
  if test "X$ac_cv_use_tpetra_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-tpetra"
  fi
else
  BUILD_TPETRA=no
fi
AM_CONDITIONAL(BUILD_TPETRA, test "X$BUILD_TPETRA" != "Xno")

if test "X$ac_cv_use_sacado" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/sacado])
  BUILD_SACADO=yes
  if test "X$ac_cv_use_sacado_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-sacado"
  fi
else
  BUILD_SACADO=no
fi
AM_CONDITIONAL(BUILD_SACADO, test "X$BUILD_SACADO" != "Xno")

if test "X$ac_cv_use_didasko" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/didasko])
  BUILD_DIDASKO=yes
  if test "X$ac_cv_use_didasko_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-didasko"
  fi
else
  BUILD_DIDASKO=no
fi
AM_CONDITIONAL(BUILD_DIDASKO, test "X$BUILD_DIDASKO" != "Xno")

if test "X$ac_cv_use_webtrilinos" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/WebTrilinos])
  BUILD_WEBTRILINOS=yes
  if test "X$ac_cv_use_webtrilinos_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-webtrilinos"
  fi
else
  BUILD_WEBTRILINOS=no
fi
AM_CONDITIONAL(BUILD_WEBTRILINOS, test "X$BUILD_WEBTRILINOS" != "Xno")

if test "X$ac_cv_use_external" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/external])
  BUILD_EXTERNAL=yes
  if test "X$ac_cv_use_external_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-external"
  fi
else
  BUILD_EXTERNAL=no
fi
AM_CONDITIONAL(BUILD_EXTERNAL, test "X$BUILD_EXTERNAL" != "Xno")

if test "X$ac_cv_use_moocho" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/moocho])
  BUILD_MOOCHO=yes
  if test "X$ac_cv_use_moocho_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-moocho"
  fi
else
  BUILD_MOOCHO=no
fi
AM_CONDITIONAL(BUILD_MOOCHO, test "X$BUILD_MOOCHO" != "Xno")

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
		Makefile.export.trilinos
		packages/Makefile])
 AC_OUTPUT
