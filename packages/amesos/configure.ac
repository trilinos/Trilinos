# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
# 
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2002) Sandia Corporation
# 
#
# Copyright (2001) Sandia Corportation. Under the terms of Contract 
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

AC_INIT(amesos, 0.1, ksstanl@sandia.gov)

# Hello World!
echo "----------------------------------------"
echo "Running Amesos Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.

AC_CONFIG_SRCDIR(/src/SuperludistOO.cpp)

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

AC_CONFIG_AUX_DIR(config)

# Configure should create src/Amesos_config.h from src/Amesos_config.h.in

AM_CONFIG_HEADER(src/Amesos_config.h:src/Amesos_config.h.in)

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

AM_INIT_AUTOMAKE(amesos, 0.1)

# Specify required version of autoconf.

AC_PREREQ(2.52)

# Specify revision number of configure.ac that was used to generate
# configure.

dnl Commented out this line for now because it causes headaches in
dnl committing.
dnl AC_REVISION($Revision$)

# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI

# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ gcc c++ cxx)
BASE_CXX=${CXX}
#AC_PROG_F77(f77 g77 f90 xlf90)
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Checks for --enable-debug, --enable-optimize, --enable-profile, --enable-purify
# Changes flags and compiler names accordingly
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_COMPILER_OPTIONS

# ------------------------------------------------------------------------
# Checks for extra flags
# ------------------------------------------------------------------------
# Check if --with-flags present, prepend any specs to FLAGS

TAC_ARG_WITH_FLAGS(ccflags, CCFLAGS)
TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
# No C or Fortran code needs to be compiled at this time
TAC_ARG_WITH_FLAGS(cflags, CFLAGS)
#TAC_ARG_WITH_FLAGS(fflags, FFLAGS)
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
# Special check for case of running on TFLOP - not needed for now
# ------------------------------------------------------------------------
#AC_MSG_CHECKING([if host is TFLOP (and therefore needs special TFLOP define statement)])
#if test "X${host}" = "Xtflop"; then
#  AC_MSG_RESULT(yes)
#  AC_DEFINE(TFLOP,,[define if compiling for ASCI Red, aka TFLOP])
#else
#  AC_MSG_RESULT(no)
#fi

# ------------------------------------------------------------------------
# Checks for header files
# ------------------------------------------------------------------------

AC_CHECK_HEADERS([cstdlib stdlib.h cstdio stdio.h cassert assert.h],,,)
AC_CHECK_HEADERS([iomanip iomanip.h iostream iostream.h string string.h],,,)

# **AMESOS SPECIFIC**
# The next two rows of headers are not verified to be present below.  At this time
# configure could complete without required headers.  We need to upgrade to include
# c std and vanilla versions of these headers where appropriate.  A centralized ConfigDefs
# file should also be written.  Errors could also occur if the c std version of a header is
# present because currently the vanilla header is simply included - no check is done yet.

AC_CHECK_HEADERS([vector map time.h fstream sys/utsname.h unistd.h],,,)
AC_CHECK_HEADERS([sys/time.h sys/resource.h],,,)

if test "X$ac_cv_header_cstdlib" != "Xyes" && test "X$ac_cv_header_stdlib_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdlib or stdlib.h.  Amesos requires at least one or the other.]])
fi

if test "X$ac_cv_header_cstdio" != "Xyes" && test "X$ac_cv_header_stdio_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdio or stdio.h.  Amesos requires at least one or the other.]])
fi

if test "X$ac_cv_header_cassert" != "Xyes" && test "X$ac_cv_header_assert_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cassert or assert.h.  Amesos requires at least one or the other.]])
fi

if test "X$ac_cv_header_string" != "Xyes" && test "X$ac_cv_header_string_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find string or string.h.  Amesos requires at least one or the other.]])
fi

if test "X$ac_cv_header_iostream" != "Xyes" && test "X$ac_cv_header_iostream_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find iostream or iostream.h.  Amesos requires at least one or the other.]])
fi

if test "X$ac_cv_header_iomanip" != "Xyes" && test "X$ac_cv_header_iomanip_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find iomanip or iomanip.h.  Amesos requires at least one or the other.]])
fi  

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE

# ------------------------------------------------------------------------
# Checks for compiler characteristics
# ------------------------------------------------------------------------
AC_CXX_NAMESPACES
AC_CXX_HAVE_STL
AC_CXX_BOOL
AC_CXX_MUTABLE
AC_CXX_NEW_FOR_SCOPING
AC_CXX_STD_SPRINTF

# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------

#TAC_ARG_ENABLE_OPTION(threads, [enable shared memory threads], THREADS, no)
#AM_CONDITIONAL(HAVE_THREADS, [test "X$ac_cv_use_threads" != "Xno"])

#TAC_ARG_ENABLE_OPTION(fatal_messages, [enable printing of fatal error messages], 
#                   FATAL_MESSAGES, yes)
#TAC_ARG_ENABLE_OPTION(warning_messages, [enable printing of warning error messages], 
#                   WARNING_MESSAGES, no)
#TAC_ARG_ENABLE_OPTION(format_io, [Use package-defined IO formatting for output], 
#                   FORMAT_IO, yes)
#TAC_ARG_WITH_PACKAGE(zoltan, [Enable Zoltan interface support], ZOLTAN, no)
#TAC_ARG_WITH_PACKAGE(rtop, [Enable RTOp interface support], RTOP, no)
#TAC_ARG_WITH_PACKAGE(esi, [Enable ESI interface support], ESI, yes)
#TAC_ARG_WITH_PACKAGE(superlu, [Enable SuperLU interface support], SUPERLU, no)

# ------------------------------------------------------------------------
# Specify other directories 
# ------------------------------------------------------------------------

# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS

# Define F77_FUNC that will be used to link with Fortran subroutines, e.g., BLAS/LAPACK
#AC_F77_WRAPPERS

# ------------------------------------------------------------------------
# Checks for addons
# ------------------------------------------------------------------------

TAC_ARG_WITH( trilinos3pldir, [Set Trilinos3PL directory], ../../../../Trilinos3PL )

## echo "  tac_with_trilinos3pldir =  ${tac_with_trilinos3pldir} " 
AC_CHECK_FILE(${tac_with_trilinos3pldir})

#
#  --enable-amesos-slud, --with-amesos-sluddir, --with-amesos-sludlib, 
#    --with-amesos-sludincdir
#
TAC_ARG_ENABLE_FEATURE(amesos_slud,[Enable Amesos interface support for SuperLU Dist], 
		AMESOS_SLUD, no)


#
#  Test to see if everything is where it ought to be:
#    Make sure that supermatrix.h is in the SuperLU serial include directory
#    Make sure that SuperLU_timer is in the SuperLU serial library
#
if test "X$ac_cv_use_amesos_slud" != "Xno"; then
  AC_LANG_SAVE
  AC_LANG_C

  TAC_ARG_WITH( amesos_sluddir, [Set SuperLU dist directory],  ${tac_with_trilinos3pldir}/SuperLU_DIST )
  AC_CHECK_FILE(${tac_with_amesos_sluddir})

  TAC_ARG_WITH( amesos_sludlib, [Set SuperLU dist library],  ${tac_with_amesos_sluddir}/libsuperludist.a )
  AC_CHECK_FILE(${tac_with_amesos_sludlib})

  TAC_ARG_WITH( amesos_sludincdir, [Set SuperLU serial include directory],  ${tac_with_amesos_sluddir}/SRC )
  AC_CHECK_FILE(${tac_with_amesos_sludincdir})

#
# This is a bad way to check for supermatrix.h, it won't pick supermatrix.h up if
# it is in the includ epath, i.e. /usr/include - see tac_arg_check_mpi for a better way
#


  AC_CHECK_FILE(${tac_with_amesos_sludincdir}/supermatrix.h)
  AC_CHECK_FILE(/home/kstanley/Trilinos3PL/SuperLU_DIST/SRC/supermatrix.h)

  AC_CHECK_FILE(/home/kstanley/Trilinos3PL/SuperLU_DIST/SRC/supermatrix.h,,
    AC_MSG_ERROR([[Cannot find supermatrix.h, specify a location using --with-amesos_sluincdir, --with-ameos_sluddir or --with-trilinos3pldir]] ))



  AC_CHECK_FILE(${tac_with_amesos_sludincdir}/supermatrix.h,,
    AC_MSG_ERROR([[Cannot find supermatrix.h in ${tac_with_amesos_sludincdir} , specify a location using --with-amesos_sluincdir, --with-ameos_sluddir or --with-trilinos3pldir]] ))


##  LDFLAGS="${tac_with_amesos_sludlib} ${LDFLAGS}"

  CXXFLAGS="-I${tac_with_amesos_sludincdir} ${CXXFLAGS}"

    SAVE_LIBS="${LIBS}"
    LIBS="${tac_with_amesos_sludlib} ${LIBS}"
    AC_TRY_LINK_FUNC(SuperLU_timer_,,
      AC_MSG_ERROR([[Cannot find SuperLU Dist library, specify a location using --with-amesos_sludlib, (ex. --with-amesos_sludlib=/home/kstanley/Trilinos3PL/SuperLU_DIST/libsuperludist.a)]]))  
    LIBS="${SAVE_LIBS}"
  AM_CONDITIONAL(HAVE_AMESOS_SLUD, true ) 
else
  AM_CONDITIONAL(HAVE_AMESOS_SLUD, false )  
fi

AC_SUBST(AMESOS_LIBS,${tac_with_amesos_sludlib}) 





## AC_ARG_ENABLE(amesos_slud,
## [AC_HELP_STRING([--enable-amesos_slud],[compile SuperLU Dist interface])],
## [ADDON_SLUD=$enableval],
## [ADDON_SLUD=yes]
## )
## AM_CONDITIONAL(BUILD_SLUDIST, test X${ADDON_SLUD} = "Xyes")

AC_ARG_ENABLE(amesos-slus,
[AC_HELP_STRING([--enable-amesos-slus],[compile SuperLU Serial interface])],
[ADDON_SLUS=$enableval],
[ADDON_SLUS=yes]
)
AM_CONDITIONAL(BUILD_SLUSERIAL, test X${ADDON_SLUS} = "Xyes")

AC_ARG_ENABLE(amesos-kundert,
[AC_HELP_STRING([--enable-amesos-kundert],[compile Kundert Sparse])],
[ADDON_KS=$enableval],
[ADDON_KS=no]
)
AM_CONDITIONAL(BUILD_KUNDERT, test X${ADDON_KS} = "Xyes")

# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

  # Find valid SuperLU Dist and Kundert libraries
### In the future we need to be able to disable dependence on Kundert lib

## if test "X$ADDON_SLUD" != "Xno"; then
##   AC_LANG_SAVE
##   AC_LANG_C
##   AC_SEARCH_LIBS(SuperLU_timer_,[superludist],,
##   AC_MSG_ERROR([[Cannot find SuperLU Dist library, specify a location using --with-ldflags (ex. --with-ldflags=-L/usr/path/lib)]]))  
## fi

dnl Replace `main' with a function in -lm:
#AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
#AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))

dnl Check for user-specified names for LAPACK and BLAS
TAC_ARG_WITH_BLASLIB
TAC_ARG_WITH_LAPACKLIB

  ACX_BLAS([], [AC_MSG_ERROR([[Cannot find blas library, specify a path using --with-ldflags=-I<DIR> (ex. --with-ldflags=-I/usr/path/lib) or a specific library using --with-blas=DIR/LIB (for example --with-blas=/usr/path/lib/libcxml.a)]])])
  ACX_LAPACK([], [AC_MSG_ERROR([[Cannot find lapack library, specify a path using --with-ldflags=-I<DIR> (ex. --with-dlflags=-I/usr/path/lib) or a specific library using --with-lapack=DIR/LIB (ex. --with-lapack=/usr/path/lib/libcxml.a)]])])

# ------------------------------------------------------------------------
# Checks for compiler characteristics
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for linker characteristics
# ------------------------------------------------------------------------

# Determine libraries needed for linking with Fortran
#AC_F77_LIBRARY_LDFLAGS

# -- DEC CXX --
dnl if test X${target_cpu} = Xalpha && test X${CXX} = Xcxx; then
dnl    XTRALDADD="cxx_repository/*.o"
dnl fi

dnl AC_MSG_CHECKING([for extra object files to be linked])
dnl if test X$XTRALDADD = X; then
dnl AC_MSG_RESULT(none)
dnl else
dnl AC_MSG_RESULT($XTRALDADD)
dnl fi

AC_SUBST(XTRALDADD)

# ------------------------------------------------------------------------
# Checks for library functions
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Setup Platform for Kundert
# ------------------------------------------------------------------------

if test "X$ADDON_KS" = "Xyes"; then

 case "$target" in
  *freebsd*) PLATFORM=freebsd
           ;;
  *linux*) PLATFORM=linux
         ;;
  *sgi*) PLATFORM=sgi
       ;;
  *powerpc-apple*) PLATFORM=ppc
      ;;
  *powerpc-ibm-aix*) PLATFORM=aix
      ;;
  *solaris*) PLATFORM=solaris
      ;;
  i386*mingw32msvc) PLATFORM=mingw
       ;;
  alpha*-dec-osf*) PLATFORM=DEC
       ;;
  *)   if test "x$PLATFORM" = "x"; then
         AC_MSG_NOTICE([Unrecognized target, using generic linux for PLATFORM])
         Xyce_ARCH=linux
       fi
       ;;
 esac

 AC_MSG_NOTICE([Adding Platform Define $PLATFORM to CCFLAGS])
 CPPFLAGS="$CPPFLAGS -D$PLATFORM"
fi

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
#		test/Makefile
AC_CONFIG_FILES([
		Makefile 
		src/Makefile
		test/Makefile
		])

AC_OUTPUT()

# Bye World!
echo "----------------------------------------"
echo "Finished Running Amesos Configure Script"
echo "----------------------------------------"
