# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Initialization 
# ------------------------------------------------------------------------



# This must be the first line in configure.ac.
# Optional 3rd argument is email address for bugs.

AC_INIT(Zoltan, 3.6, lriesen@sandia.gov)

echo "----------------------------------------"
echo "Running Zoltan Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.  Any file in that directory will do,
# preferably one that is unlikely to be removed or renamed.

AC_CONFIG_SRCDIR([src/include/zoltan.h])

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

AC_CONFIG_AUX_DIR(config)

#
# We don't want people to configure in the source directory.  Some
# things may break.
#

if test -e configure.ac ; then
  echo "You are trying to run configure in the source directory.  This is not allowed.  Please run configure from a separate build directory."
  exit
fi

# Configure should create src/Zoltan_config.h from src/Zoltan_config.h.in
AM_CONFIG_HEADER(src/include/Zoltan_config.h:src/include/Zoltan_config.h.in)

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

# Note in header file and Makefile conditional what the host OS is
  
AM_CONDITIONAL(HOST_CYGWIN, false)
AM_CONDITIONAL(HOST_LINUX, false)
AM_CONDITIONAL(HOST_SOLARIS, false)

case $host_os in
  cygwin)
    AM_CONDITIONAL(HOST_CYGWIN, true)
    AC_DEFINE(HOST_CYGWIN,1,[software host will be cygwin])
    ;;
  linux*)
    AM_CONDITIONAL(HOST_LINUX, true)
    AC_DEFINE(HOST_LINUX,1,[software host will be linux])
    ;;
  solaris*)
    AM_CONDITIONAL(HOST_SOLARIS, true)
    AC_DEFINE(HOST_SOLARIS,1,[software host will be solaris])
    ;;
esac

# Use automake

#  - Required version of automake.
AM_INIT_AUTOMAKE(1.9.6 no-define tar-ustar)

# Specify required version of autoconf.

AC_PREREQ(2.59)

#TAC_ARG_ENABLE_OPTION(fortran, [enable Fortran support], FORTRAN_SUPPORT, no)
#This option is not currently available
TAC_ARG_ENABLE_OPTION(f90interface, [enable Fortran 90 interface (automatically enables Fortran support)], F90INTERFACE, no)
AM_CONDITIONAL(BUILD_ZOLTAN_F90_INTERFACE, [test "X$ac_cv_use_f90interface" != "Xno"])
#AM_CONDITIONAL(USE_FORTRAN, [test "X$ac_cv_use_fortran" != "Xno"])

if test "X$ac_cv_use_f90interface" = "Xyes"; then
  ac_cv_use_fortran=yes
else
  ac_cv_use_fortran=no
fi

if test "X$ac_cv_use_fortran" = "Xyes"; then
  ac_cv_use_fortran90=yes
#  AX_F90_MODULE_FLAG
#  if test "X$ax_cv_f90_modflag" = "Xunknown" ; then
#    AC_MSG_ERROR([unable to find f90 modules extension])
#  else
#    FCFLAGS="$ax_cv_f90_modflag../ $ax_cv_f90_modflag. ${FCFLAGS}"
#  fi
else
  ac_cv_use_fortran90=no
fi

TAC_ARG_ENABLE_FEATURE_SUB(zoltan, cppdriver, [Enable Zoltan's C++ driver], ZOLTAN_CPPDRIVER, yes)
AM_CONDITIONAL(HAVE_ZOLTAN_CPPDRIVER, [test "X$ac_cv_use_zoltan_cppdriver" != "Xno"])

#This can be removed after we retire the old build system
#AC_DEFINE([TRILINOS_CONFIG_H],,[Define when using the autotools to build Zoltan])

# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

# We may want to handle this differently because Zoltan requires MPI
#TAC_ARG_CONFIG_MPI

ZAC_ARG_CONFIG_MPI

#  #np# - can eliminate compiler checks below if your package does not use the
#         language corresponding to the check.  Please note that if you use
#	  F77_FUNC to determine Fortran name mangling, you should not remove
#	  the Fortran compiler check or the check for Fortran flags.  Doing
#	  so will prevent the detection of the proper name mangling in some
#	  cases.

AC_ARG_ENABLE(mpi-recv-limit,
 AS_HELP_STRING([--enable-mpi-recv-limit],[Set to the limit on the number of simultaneous MPI posted receives, if any; default is --enable-mpi-recv-limit=10]),
 [ AC_MSG_NOTICE(Building Zoltan to observe a limit of $enable_mpi_recv_limit simultaneously posted MPI receives)
   AC_SUBST(MPI_RECV_LIMIT_FLAG, "-DMPI_RECV_LIMIT=$enable_mpi_recv_limit")
 ]
 )

# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

AM_PROG_CC_C_O(cc gcc)
#Conditional only?
#if test "X$ac_cv_use_zoltan_cppdriver" = "Xyes"; then
AC_PROG_CXX(CC g++ c++ cxx)
#fi

if test "X$ac_cv_use_fortran90" = "Xyes"; then
AC_PROG_FC(FC gfortran f90 xlf90 f95)
fi

if test "X$ac_cv_use_fortran" = "Xyes"; then
AC_PROG_FC(FC gfortran f90 xlf90 f95)
fi

AC_PROG_RANLIB

# Check if --with-flags present, prepend any specs to FLAGS

TAC_ARG_WITH_FLAGS(ccflags, CCFLAGS)
#if test "X$ac_cv_use_zoltan_cppdriver" = "Xyes"; then
TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
CXXFLAGS="-DMPICH_IGNORE_CXX_SEEK  ${CXXFLAGS}"
#fi
TAC_ARG_WITH_FLAGS(cflags, CFLAGS)
TAC_ARG_WITH_FLAGS(fcflags, FCFLAGS)
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
# Add this later
TAC_ARG_ENABLE_EXPORT_MAKEFILES(yes)

# ------------------------------------------------------------------------
# Checks for special package flags 
# ------------------------------------------------------------------------

ZAC_ARG_WITH_ID()

TAC_ARG_ENABLE_OPTION(gzip, [enable zlib support for driver], GZIP, no)
if test "X$ac_cv_use_gzip" != "Xno"; then
   LIBS="-lz ${LIBS}"
fi
AM_CONDITIONAL(BUILD_GZIP, [test "X$ac_cv_use_gzip" != "Xno"])

TAC_ARG_WITH_PACKAGE(parmetis, [Enable Parmetis support.], PARMETIS, no)
AM_CONDITIONAL(BUILD_PARMETIS, [test "X$ac_cv_use_parmetis" != "Xno"])
TAC_ARG_WITH_3PL_SUB(parmetis, libdir, [Specify where the Parmetis library is located.  Ex. /path/to/library])
if test "X$tac_with_parmetis_libdir" != "Xno"; then
  LIBS="-L${tac_with_parmetis_libdir} -lparmetis -lmetis ${LIBS}"
fi

TAC_ARG_WITH_3PL_SUB(parmetis, incdir, [Specify where the Parmetis header files are located.  Ex. /path/to/headers])
# It was necessary to move the parmetis include directories before the scotch
# include directories.
#if test "X$tac_with_parmetis_incdir" != "Xno"; then
#  CPPFLAGS="-I${tac_with_parmetis_incdir}  ${CPPFLAGS}"
#fi


TAC_ARG_WITH_PACKAGE(scotch, [Enable Scotch support.], SCOTCH, no)
AM_CONDITIONAL(BUILD_SCOTCH, [test "X$ac_cv_use_scotch" != "Xno"])
TAC_ARG_WITH_3PL_SUB(scotch, libdir, [Specify where the Scotch library is located.  Ex. /path/to/library])
if test "X$tac_with_scotch_libdir" != "Xno"; then
  LIBS="-L${tac_with_scotch_libdir} -lptscotch  -lptscotcherr -lptscotcherrexit -lscotch -lscotcherr -lscotcherrexit ${LIBS}"
fi
TAC_ARG_WITH_3PL_SUB(scotch, incdir, [Specify where the Scotch header files are located.  Ex. /path/to/headers])
if test "X$tac_with_scotch_incdir" != "Xno"; then
  CPPFLAGS="-I${tac_with_scotch_incdir}  ${CPPFLAGS}"
fi
# The parmetis include directories need to be before the Scotch include
# directories.
if test "X$tac_with_parmetis_incdir" != "Xno"; then
  CPPFLAGS="-I${tac_with_parmetis_incdir}  ${CPPFLAGS}"
fi

TAC_ARG_WITH_PACKAGE(patoh, [Enable Patoh support.], PATOH, no)
TAC_ARG_WITH_3PL_SUB(patoh, libdir, [Specify where the Patoh library is located.  Ex. /path/to/library])
if test "X$tac_with_patoh_libdir" != "Xno"; then
  LIBS="-L${tac_with_patoh_libdir} -lpatoh ${LIBS}"
fi
TAC_ARG_WITH_3PL_SUB(patoh, incdir, [Specify where the Patoh header files are located.  Ex. /path/to/headers])
if test "X$tac_with_patoh_incdir" != "Xno"; then
  CPPFLAGS="-I${tac_with_patoh_incdir}  ${CPPFLAGS}"
fi

#TAC_ARG_WITH_PACKAGE(drum, [Enable Drum support.  Library name as well as library and include paths must be specified using LDFLAGS and CPPFLAGS.], DRUM, no)
#TAC_ARG_WITH_3PL_SUB(drum, libdir, [Specify where the Drum library is located.  Ex. /path/to/library])
#if test "X$tac_with_drum_libdir" != "Xno"; then
#  LIBS="-L${tac_with_drum_libdir} ${LIBS}"
#fi
#TAC_ARG_WITH_3PL_SUB(drum, incdir, [Specify where the Drum header files are located.  Ex. /path/to/headers])
#if test "X$tac_with_drum_incdir" != "Xno"; then
#  CPPFLAGS="-I${tac_with_drum_incdir}  ${CPPFLAGS}"
#fi

TAC_ARG_WITH_PACKAGE(nemesis_exodus, [Enable Nemesis/Exodus support for the Zdrive test executable.  Library and include paths must be specified using LDFLAGS and CFLAGS.], NEMESIS_EXODUS, no)
#I am not sure if we can do this.  We may not be able to get things in the right
#order. JW.
if test "X$ac_cv_use_nemesis_exodus" != "Xno"; then
  LIBS="-lnemIc -lexoIIv2c -lnetcdf ${LIBS}"
fi

LIBS="${LIBS} -lm"

#TAC_ARG_WITH_PACKAGE(parkway, [Enable Parkway support.], PARKWAY, no)
#TAC_ARG_WITH_3PL_SUB(parkway, libdir, [Specify where the Parkway library is located.  Ex. /path/to/library])
#if test "X$tac_with_parkway_libdir" != "Xno"; then
#  LIBS="-L${tac_with_parkway_libdir} -lparkway -lpmpich++ -lstdc++ ${LIBS}"
#fi
#TAC_ARG_WITH_3PL_SUB(parkway, incdir, [Specify where the Parkway header files are located.  Ex. /path/to/headers])
#if test "X$tac_with_parkway_incdir" != "Xno"; then
#  CPPFLAGS="-I${tac_with_parkway_incdir}  ${CPPFLAGS}"
#fi

# ------------------------------------------------------------------------
# Checks if tests and examples should be built
# ------------------------------------------------------------------------


TAC_ARG_ENABLE_FEATURE(tests, [Make tests for all Trilinos packages buildable with 'make tests'], TESTS, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( zoltan, tests, [Make Zoltan tests buildable with 'make tests'], ZOLTAN_TESTS)
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_zoltan_tests" != "Xno")

TAC_ARG_ENABLE_FEATURE(examples, [Make examples for all Trilinos packages buildable with 'make examples'], EXAMPLES, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( zoltan, examples, [Make Zoltan examples buildable with 'make examples'], ZOLTAN_EXAMPLES)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_zoltan_examples" != "Xno")

#We now build tests and examples through separate make targets, rather than
#during "make".  We still need to conditionally include the test and example
#in SUBDIRS, even though SUB_TEST and SUB_EXAMPLE will never be
#defined, so that the tests and examples are included in the distribution
#tarball.
AM_CONDITIONAL(SUB_TEST, test "X$ac_cv_use_sub_test" = "Xyes")
AM_CONDITIONAL(SUB_EXAMPLE, test "X$ac_cv_use_sub_example" = "Xyes")

#TAC_ARG_ENABLE_FEATURE(libcheck, [Check for some third-party libraries.  (Cannot be disabled unless tests and examples are also disabled.)], LIBCHECK, yes)

# ------------------------------------------------------------------------
# Specify other directories 
# ------------------------------------------------------------------------

# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Iincdir1 -Iincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS

# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

ax_cv_f90_modulecase="lower"
FC_VENDOR="none"

# Define F77_FUNC that will be used to link with Fortran subroutines.
if test "X$ac_cv_use_fortran" != "Xno"; then
  AC_FC_WRAPPERS
  AC_FC_SRCEXT(f)
  AC_FC_SRCEXT(f90)
  # Determine compile-line flag for F90 modules (e.g., -M).
  AX_F90_MODULE_FLAG
  if test "X$ax_cv_f90_modflag" = "Xunknown" ; then
    AC_MSG_ERROR([unable to find f90 modules extension])
  else
    FCFLAGS="$ax_cv_f90_modflag../ $ax_cv_f90_modflag. ${FCFLAGS}"
  fi
  # Determine case (upper or lower) of F90 module files. 
  AX_F90_MODULE_CASE
  echo "KDDKDD " $ax_cv_f90_modulecase 

  WK_FC_GET_VENDOR()
  echo "KDDKDD " $FC_VENDOR
  if test "$FC_VENDOR" = "PGI"; then
    FCFLAGS="-DPGI ${FCFLAGS}"
    CFLAGS="-DPGI ${CFLAGS}"
  fi
  if test "$FC_VENDOR" = "Lahey"; then
    FCFLAGS="-DFUJITSU ${FCFLAGS}"
    CFLAGS="-DFUJITSU ${CFLAGS}"
  fi
fi
AM_CONDITIONAL(F90_MODULE_LOWERCASE, [test "X$ax_cv_f90_modulecase" != "Xupper"])
AM_CONDITIONAL(NAG_F90_COMPILER, [test "X$FC_VENDOR" == "XNAG"])

# If tests, examples and libcheck are disabled, we don't have to check
# for these libraries.

#if test "X$ac_cv_use_new_package_examples" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then
#if test "X$ac_cv_use_zoltan_tests" != "Xno" || test "X$ac_cv_use_zoltan_examples" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then
# checks for libraries now
dnl Replace `main' with a function in -lm:
#AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
#AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))
#fi
# end of the list of libraries that don't need to be checked for if
# tests and examples are disabled.

# ------------------------------------------------------------------------
# Checks for linker characteristics
# ------------------------------------------------------------------------

# Determine libraries needed for linking with Fortran
#AC_FC_LIBRARY_LDFLAGS
if test "X$ac_cv_use_fortran" = "Xyes"; then
AC_FC_LIBRARY_LDFLAGS
fi

# ------------------------------------------------------------------------
# Zoltan may be built via autotools, cmake, or our native makefile
# ------------------------------------------------------------------------

CPPFLAGS="-DAUTOTOOLS_BUILD  ${CPPFLAGS}"

# ------------------------------------------------------------------------
# Perform substitutions in output files
# ------------------------------------------------------------------------

AC_SUBST(ac_aux_dir)

# ------------------------------------------------------------------------
# Output files
# ------------------------------------------------------------------------
#
AC_CONFIG_FILES([
		Makefile
		src/Makefile
		src/driver/Makefile
		src/fdriver/Makefile
		Makefile.export.zoltan
		example/Makefile
		example/C/Makefile
		example/CPP/Makefile
		siMPI/Makefile
		siMPI/pyMPI/Makefile
		siMPI/pyMPI/siMPI/Makefile
		])
#		test/Large_Data/Makefile
#		example/lib/Makefile
#		example/C/Makefile
#		example/C/sparse_matrix/Makefile
#		example/CPP/Makefile
#		Utilities/Makefile

AC_OUTPUT()

echo "---------------------------------------------"
echo "Finished Running Zoltan Configure Script"
echo "---------------------------------------------"

