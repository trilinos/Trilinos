# @HEADER
# @HEADER

#
#  Refer to packages/new_package/README for suggestions on 
#  adding a new package.  
#

# Process this file with autoconf to produce a configure script.
AC_INIT(Trilinos, 10.0d, maherou@sandia.gov)
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

AM_INIT_AUTOMAKE(1.10 no-define tar-ustar)

# Specify required version of autoconf.

AC_PREREQ(2.61)

# Must indicate for packages that can be configured from the package level that the 
# configure/build process is happening within Trilinos.

ac_configure_args="$ac_configure_args --with-trilinos"

echo "NOTE: If a package is not enabled based on the"
echo "configure options, but is required by another"
echo "package that has been enabled, it will"
echo "automatically be configured and built, even if"
echo "the list below does not indicate so."

TAC_ARG_ENABLE_FEATURE_SUB( default, packages, [Build all default packages unless explicitly disabled.  When this option is disabled, only packages explicitly enabled and packages needed to satisfy the dependencies of enabled packages will be built], DEFAULT_PACKAGES, yes)
 
TAC_ARG_ENABLE_FEATURE(amesos, [Build amesos, a direct sparse solver package], AMESOS, no)

TAC_ARG_ENABLE_FEATURE(anasazi, [Build anasazi, an eigenvalue solver package], ANASAZI, no)

TAC_ARG_ENABLE_FEATURE(aristos, [Build aristos, an optimization package], ARISTOS, no)

TAC_ARG_ENABLE_DEFAULT_FEATURE(aztecoo, [Build aztecoo, a linear solver and preconditioner package], AZTECOO)

TAC_ARG_ENABLE_FEATURE(belos, [Build belos, a linear solver package], BELOS, no)

TAC_ARG_ENABLE_FEATURE(claps, [Build claps, a domain decomposition preconditioner package], CLAPS, no)

TAC_ARG_ENABLE_FEATURE(didasko, [Build didasko, a Trilinos tutorial package], DIDASKO, no)

TAC_ARG_ENABLE_DEFAULT_FEATURE(epetra, [Build epetra, a basic linear algebra package], EPETRA)

TAC_ARG_ENABLE_DEFAULT_FEATURE(epetraext, [Build epetraext, an epetra extensions package], EPETRAEXT)

#The external package is a mechanism for adding external packages to the Trilinos
#build process.
TAC_ARG_ENABLE_FEATURE(external, [Build external package in packages/external], EXTERNAL, no)

TAC_ARG_ENABLE_DEFAULT_FEATURE(ifpack, [Build ifpack, an algebraic preconditioner package], IFPACK)

TAC_ARG_ENABLE_FEATURE(fei, [Build fei, a sparse linear-system assembly package], FEI, no)

TAC_ARG_ENABLE_FEATURE(ctrilinos, [Build CTrilinos, a package providing C wrappers for other Trilinos packages], CTRILINOS, no)

TAC_ARG_ENABLE_FEATURE(fortrilinos, [Build ForTrilinos, a package providing Fortran wrappers for other Trilinos packages], FORTRILINOS, no)

TAC_ARG_ENABLE_FEATURE(intrepid, [Build intrepid, a discretization package], INTREPID, no)

TAC_ARG_ENABLE_FEATURE(isorropia, [Build isorropia, a repartitioning/rebalancing package], ISORROPIA, no)

TAC_ARG_ENABLE_FEATURE(kokkos, [Build kokkos, a sparse kernels package], KOKKOS, no)

TAC_ARG_ENABLE_FEATURE(komplex, [Build komplex, a complex linear solver package], KOMPLEX, no)

TAC_ARG_ENABLE_DEFAULT_FEATURE(loca, [Build loca, a bifurcation analysis package], LOCA)

TAC_ARG_ENABLE_FEATURE(meros, [Build meros, a segregated preconditioner package], MEROS, no)

TAC_ARG_ENABLE_FEATURE(moertel, [Build moertel, a mortar methods package], MOERTEL, no)

TAC_ARG_ENABLE_DEFAULT_FEATURE(ml, [Build ml, a multi-level preconditioner package], ML)

## #np# - Copy the next line, and change the name.  Place the line in 
## #np# - alphabetical order based on package name.
TAC_ARG_ENABLE_DEFAULT_FEATURE(new_package, [Build new_package, a working package prototype], NEW_PACKAGE)

TAC_ARG_ENABLE_DEFAULT_FEATURE(nox, [Build nox, a nonlinear solver package], NOX)

TAC_ARG_ENABLE_FEATURE(phdmesh, [Build phdMesh, a Parallel Heterogneous Dynamic unstructured Mesh package], PHDMESH, no)

TAC_ARG_ENABLE_FEATURE(pliris, [Build pliris, a dense direct solver package], PLIRIS, no)

TAC_ARG_ENABLE_FEATURE(pytrilinos, [Build PyTrilinos, a Python wrapper package], PYTRILINOS, no)

TAC_ARG_ENABLE_FEATURE(rbgen, [Build rbgen], RBGEN, no)

TAC_ARG_ENABLE_FEATURE(rythmos, [Build rythmos, a transient integration package], RYTHMOS, no)

TAC_ARG_ENABLE_FEATURE(sacado, [Build sacado, an automatic differentation package], SACADO, no)

TAC_ARG_ENABLE_FEATURE(shards, [Build shards, a package of shared discretization tools], SHARDS, no)

TAC_ARG_ENABLE_FEATURE(stokhos, [Build stokhos, an intrusive stochastic Galerkin package], STOKHOS, no)

TAC_ARG_ENABLE_FEATURE(sundance, [Build Sundance, a system for development of parallel finite-element solutions of PDEs], SUNDANCE, no)

TAC_ARG_ENABLE_FEATURE(teuchos, [Build teuchos, a tools package], TEUCHOS, no)

TAC_ARG_ENABLE_FEATURE(pamgen, [Build pamgen, a tools package], PAMGEN, no)

TAC_ARG_ENABLE_FEATURE(rtop, [Build RTOp, a vector reduction/transformation operator package], RTOP, no)

TAC_ARG_ENABLE_FEATURE(threadpool, [Build ThreadPool], THREADPOOL, no)

TAC_ARG_ENABLE_FEATURE(thyra, [Build thyra, an abstract linear solver package], THYRA, no)

TAC_ARG_ENABLE_FEATURE(stratimikos, [Build stratimikos, a linear solver strategy package], STRATIMIKOS, no)

TAC_ARG_ENABLE_FEATURE(moocho, [Build moocho, an optimization package], MOOCHO, no)

TAC_ARG_ENABLE_FEATURE(tpetra, [Build tpetra, a templated linear algebra package], TPETRA, no)

TAC_ARG_ENABLE_FEATURE(trikota, [Build TriKota, a package providing Trilinos/Dakota interoperability], TRIKOTA, no)

TAC_ARG_ENABLE_FEATURE(trilinoscouplings, [Build trilinoscouplings, a package of interfaces between Trilinos packages], TRILINOSCOUPLINGS, no)

TAC_ARG_ENABLE_FEATURE(triutils, [Build triutils, a utility package], TRIUTILS, no)

TAC_ARG_ENABLE_FEATURE(galeri, [Build galeri, a gallery package], GALERI, no)

TAC_ARG_ENABLE_FEATURE(webtrilinos, [Build WebTrilinos, a web interface package], WEBTRILINOS, no)

TAC_ARG_ENABLE_FEATURE(zoltan, [Build zoltan, a partitioning and load balancing package], ZOLTAN, no)

TAC_ARG_ENABLE_FEATURE(shared, [Build shared libraries in addition to static libraries], SHARED, no)

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
 	test "X$ac_cv_use_aristos" != "Xno" ||
	test "X$ac_cv_use_intrepid" != "Xno" ||
	test "X$ac_cv_use_aztecoo_teuchos" != "Xno" ||
	test "X$ac_cv_use_epetraext" != "Xno" ||
	test "X$ac_cv_use_meros" != "Xno" ||
	test "X$ac_cv_use_moocho" != "Xno" ||
	test "X$ac_cv_use_nox" != "Xno" ||
	test "X$ac_cv_use_rbgen" != "Xno" ||
	test "X$ac_cv_use_rtop" != "Xno" ||
	test "X$ac_cv_use_rythmos" != "Xno" ||
	test "X$ac_cv_use_sacado" != "Xno" ||
	test "X$ac_cv_use_stokhos" != "Xno" ||
	test "X$ac_cv_use_sundance" != "Xno" ||
	test "X$ac_cv_use_isorropia" != "Xno" ||
	test "X$ac_cv_use_moertel" != "Xno" ||
	test "X$ac_cv_use_webtrilinos" != "Xno" ||
	test "X$ac_cv_use_ifpack" != "Xno" ||
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


if test "X$ac_cv_use_pamgen" != "Xno" ; then
  AC_CONFIG_SUBDIRS([packages/pamgen])
  BUILD_PAMGEN=yes
  if test "X$ac_cv_use_pamgen_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-teuchos"
  fi
else
  BUILD_PAMGEN=no
fi
AM_CONDITIONAL(BUILD_PAMGEN, test "X$BUILD_PAMGEN" != "Xno")





if test "X$ac_cv_use_rtop" != "Xno" ||
	test "X$ac_cv_use_moocho" != "Xno" ||
	test "X$ac_cv_use_meros" != "Xno" ||
	test "X$ac_cv_use_rythmos" != "Xno" ||
	test "X$ac_cv_use_stratimikos" != "Xno" ||
	test "X$ac_cv_use_sundance" != "Xno" ||
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

if test "X$ac_cv_use_threadpool" != "Xno" ||
	test "X$ac_cv_use_phdmesh" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/ThreadPool])
  BUILD_THREADPOOL=yes
  if test "X$ac_cv_use_threadpool_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-threadpool"
  fi
else
  BUILD_THREADPOOL=no
fi
AM_CONDITIONAL(BUILD_THREADPOOL, test "X$BUILD_THREADPOOL" != "Xno")

if test "X$ac_cv_use_thyra" != "Xno" ||
	test "X$ac_cv_use_meros" != "Xno" ||
	test "X$ac_cv_use_rythmos" != "Xno" ||
	test "X$ac_cv_use_sundance" != "Xno" ||
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

if test "X$ac_cv_use_stratimikos" != "Xno" ||
	test "X$ac_cv_use_sundance" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/stratimikos])
  BUILD_STRATIMIKOS=yes
  if test "X$ac_cv_use_stratimikos_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-stratimikos"
  fi
else
  BUILD_STRATIMIKOS=no
fi
AM_CONDITIONAL(BUILD_STRATIMIKOS, test "X$BUILD_STRATIMIKOS" != "Xno")

if test "X$ac_cv_use_kokkos" != "Xno"; then
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
	test "X$ac_cv_use_aristos" != "Xno" ||
        test "X$ac_cv_use_claps" != "Xno" ||
        test "X$ac_cv_use_isorropia" != "Xno" ||
	test "X$ac_cv_use_meros" != "Xno" ||
        test "X$ac_cv_use_moertel" != "Xno" ||
	test "X$ac_cv_use_sundance" != "Xno" ||
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

if test "X$ac_cv_use_fei" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/fei])
  BUILD_FEI=yes
else
  BUILD_FEI=no
fi
AM_CONDITIONAL(BUILD_FEI, test "X$BUILD_FEI" != "Xno")

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
	test "X$ac_cv_use_sundance" != "Xno" ||
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
	test "X$ac_cv_use_sundance" != "Xno" ||
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
	test "X$ac_cv_use_sundance" != "Xno" ||
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
	test "X$ac_cv_use_sundance" != "Xno" ||
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
	test "X$ac_cv_use_sundance" != "Xno" ||
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

if test "X$ac_cv_use_rbgen" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/rbgen])
  BUILD_RBGEN=yes
  if test "X$ac_cv_use_rbgen_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-rbgen"
  fi
else
  BUILD_RBGEN=no
fi
AM_CONDITIONAL(BUILD_RBGEN, test "X$BUILD_RBGEN" != "Xno")

if test "X$ac_cv_use_aristos" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/aristos])
  BUILD_ARISTOS=yes
  if test "X$ac_cv_use_aristos_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-aristos"
  fi
else
  BUILD_ARISTOS=no
fi
AM_CONDITIONAL(BUILD_ARISTOS, test "X$BUILD_ARISTOS" != "Xno")

if test "X$ac_cv_use_intrepid" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/intrepid])
  BUILD_INTREPID=yes
  if test "X$ac_cv_use_intrepid_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-intrepid"
  fi
else
  BUILD_INTREPID=no
fi
AM_CONDITIONAL(BUILD_INTREPID, test "X$BUILD_INTREPID" != "Xno")

if test "X$ac_cv_use_shards" != "Xno" ||
	test "X$ac_cv_use_intrepid" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/shards])
  BUILD_SHARDS=yes
  if test "X$ac_cv_use_shards_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-shards"
  fi
else
  BUILD_SHARDS=no
fi
AM_CONDITIONAL(BUILD_SHARDS, test "X$BUILD_SHARDS" != "Xno")

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
	test "X$ac_cv_use_sundance" != "Xno" ||
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

if test "X$ac_cv_use_stokhos" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/stokhos])
  BUILD_STOKHOS=yes
  if test "X$ac_cv_use_stokhos_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-stokhos"
  fi
else
  BUILD_STOKHOS=no
fi
AM_CONDITIONAL(BUILD_STOKHOS, test "X$BUILD_STOKHOS" != "Xno")

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

if test "X$ac_cv_use_nox" != "Xno" ||
	test "X$ac_cv_use_sundance" != "Xno" ||
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

if test "X$ac_cv_use_phdmesh" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/phdmesh])
  BUILD_PHDMESH=yes
  if test "X$ac_cv_use_phdmesh_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-phdmesh"
  fi
else
  BUILD_PHDMESH=no
fi
AM_CONDITIONAL(BUILD_PHDMESH, test "X$BUILD_PHDMESH" != "Xno")

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

if test "X$ac_cv_use_sundance" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/Sundance])
  BUILD_SUNDANCE=yes
  if test "X$ac_cv_use_sundance_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-sundance"
  fi
else
  BUILD_SUNDANCE=no
fi
AM_CONDITIONAL(BUILD_SUNDANCE, test "X$BUILD_SUNDANCE" != "Xno")

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
	test "X$ac_cv_use_shared" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/PyTrilinos])
  BUILD_PYTRILINOS=yes
#Below we have to list all packages that depend on pytrilinos
#because if pytrilinos is being configured only because shared is
#enabled, we don't need to enable pytrilinos.  Note that adding
#additional dependent packages would result in a statement like
#if pytrilinos_explicit != yes && (webtrilinos != no || packageB != no)
  if test "X$ac_cv_use_pytrilinos_explicit" != "Xyes" &&
	test "X$ac_cv_use_webtrilinos" != "Xno"; then
    ac_configure_args="$ac_configure_args --enable-pytrilinos"
  fi
  if test "X$ac_cv_use_shared" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-shared"
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

if test "X$ac_cv_use_trikota" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/TriKota])
  BUILD_TRIKOTA=yes
  if test "X$ac_cv_use_trikota_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-trikota"
  fi
else
  BUILD_TRIKOTA=no
fi
AM_CONDITIONAL(BUILD_TRIKOTA, test "X$BUILD_TRIKOTA" != "Xno")

if test "X$ac_cv_use_trilinoscouplings" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/trilinoscouplings])
  BUILD_TRILINOSCOUPLINGS=yes
  if test "X$ac_cv_use_trilinoscouplings_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-trilinoscouplings"
  fi
else
  BUILD_TRILINOSCOUPLINGS=no
fi
AM_CONDITIONAL(BUILD_TRILINOSCOUPLINGS, test "X$BUILD_TRILINOSCOUPLINGS" != "Xno")

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

if test "X$ac_cv_use_ctrilinos" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/CTrilinos])
  BUILD_CTRILINOS=yes
  if test "X$ac_cv_use_ctrilinos_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-ctrilinos"
  fi
else
  BUILD_CTRILINOS=no
fi
AM_CONDITIONAL(BUILD_CTRILINOS, test "X$BUILD_CTRILINOS" != "Xno")

if test "X$ac_cv_use_zoltan" != "Xno" ||
        test "X$ac_cv_use_isorropia" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/zoltan])
  BUILD_ZOLTAN=yes
  if test "X$ac_cv_use_zoltan_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-zoltan"
  fi
else
  BUILD_ZOLTAN=no
fi
AM_CONDITIONAL(BUILD_ZOLTAN, test "X$BUILD_ZOLTAN" != "Xno")

if test "X$ac_cv_use_fortrilinos" != "Xno"; then
  AC_CONFIG_SUBDIRS([packages/ForTrilinos])
  BUILD_FORTRILINOS=yes
  if test "X$ac_cv_use_fortrilinos_explicit" != "Xyes"; then
    ac_configure_args="$ac_configure_args --enable-fortrilinos"
  fi
else
  BUILD_FORTRILINOS=no
fi
AM_CONDITIONAL(BUILD_FORTRILINOS, test "X$BUILD_FORTRILINOS" != "Xno")

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

The options listed above are specific to the top Trilinos
level.  These options are limited to deciding which packages
should be built.  All other options such as compilers and
compiler flags are set at the package level.  For a list of
all options for all packages, type "configure --help=recursive".
The output for this command is bulky.  Consider using the
"more" or "grep" command to find a specific option.

To view a list of options for a specific package, cd to 
"Trilinos/package/package_name" (replacing "package_name"
with the name of the package) and type "configure --help".

THE TRILINOS USER GUIDE contains a lot of additional information
about configuring and building Trilinos, as well as other
topics that will be helpful for Trilinos Users.  The Trilinos
User Guide, as well as other useful documents, can be found at:

  http://trilinos.sandia.gov/documentation.html

])
AC_CONFIG_FILES([
		Makefile
		Makefile.export.trilinos
		packages/Makefile])
 AC_OUTPUT
