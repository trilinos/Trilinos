# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
# 
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2002,2003) Sandia Corporation
# 
#
# Copyright (2001,2002,2003) Sandia Corportation. Under the terms of Contract 
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

# package name, version number, and e-mail address below
AC_INIT(ml, 3.0, jhu@sandia.gov)

# Hello World!
#
echo "----------------------------------------"
echo "Running ML Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.  Any file in that directory will do,
# preferably one that is unlikely to be removed or renamed.

AC_CONFIG_SRCDIR(/src/Main/ml_struct.c)

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

#
# This a a subdirectory of Trilinos/packages/ml that contains
# a number of useful m4 files.  Any package that uses these m4 files
# will have to add them to the repository.  The m4 files perform various
# common checks and set common variable values to minimize redundant work.
#
AC_CONFIG_AUX_DIR(config)
# Configure should create src/ml_config.h from src/ml_config.h.in

AM_CONFIG_HEADER(src/ml_config.h:src/ml_config.h.in)

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

# Package name and version below
AM_INIT_AUTOMAKE(ml, 3.0)

# Specify required version of autoconf.

AC_PREREQ(2.52)

#  can eliminate compiler checks below if your package does not use the
#  language corresponding to the check
# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI

AC_PROG_CXX(CC g++ gcc c++ cxx)

TAC_ARG_ENABLE_OPTION(ml_multiple_rhs, [experimental], ML_MULTIPLE_RHS, no)
#TAC_ARG_ENABLE_OPTION(ml_multiple_rhs, [Enable optimized code for multiple right hand sides.   Note: this will automatically cause all ML object files to be built with the c++ compiler.], ML_MULTIPLE_RHS, no)
if test "X$ac_cv_use_ml_multiple_rhs" != "Xno"; then
CC=${CXX}
fi
AC_PROG_CC(cc gcc)
AC_PROG_F77(f77 g77 f90 xlf90)
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------



# ------------------------------------------------------------------------
# Checks for extra flags
# You can pull flag checks for languages not used in package 
# ------------------------------------------------------------------------
# Check if --with-flags present, prepend any specs to FLAGS

TAC_ARG_WITH_FLAGS(ccflags, CCFLAGS)
TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
TAC_ARG_WITH_FLAGS(cflags, CFLAGS)
TAC_ARG_WITH_FLAGS(fflags, FFLAGS)
TAC_ARG_WITH_LIBS(libs, LIBS)
TAC_ARG_WITH_FLAGS(ldflags, LDFLAGS)

# ------------------------------------------------------------------------
# Alternate archiver
# ------------------------------------------------------------------------

TAC_ARG_WITH_AR

# ------------------------------------------------------------------------
# MPI link check
# ------------------------------------------------------------------------
TAC_ARG_CHECK_MPI

# Customize the checks for headers needed in package being added.
# NOTE: This section still needs to be updated for ml, it is still 
# code from Epetra
# ------------------------------------------------------------------------
# Checks for header files - used in Epetra_ConfigDefs.h 
# ------------------------------------------------------------------------

AC_CHECK_HEADERS([cstdlib stdlib.h cstdio stdio.h cassert assert.h],,,)
AC_CHECK_HEADERS([cmath math.h iomanip iomanip.h iostream iostream.h],,,)
AC_CHECK_HEADERS([string string.h],,,)

if test "X$ac_cv_header_cstdlib" != "Xyes" && test "X$ac_cv_header_stdlib_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdlib or stdlib.h.  ml requires at least one or the other.]])
fi

if test "X$ac_cv_header_cstdio" != "Xyes" && test "X$ac_cv_header_stdio_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdio or stdio.h.  ml requires at least one or the other.]])
fi

if test "X$ac_cv_header_cassert" != "Xyes" && test "X$ac_cv_header_assert_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cassert or assert.h.  ml requires at least one or the other.]])
fi

if test "X$ac_cv_header_string" != "Xyes" && test "X$ac_cv_header_string_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find string or string.h.  ml requires at least one or the other.]])
fi

if test "X$ac_cv_header_iostream" != "Xyes" && test "X$ac_cv_header_iostream_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find iostream or iostream.h.  ml requires at least one or the other.]])
fi

if test "X${host}" != "Xtflop"; then
  if test "X$ac_cv_header_cmath" != "Xyes" && test "X$ac_cv_header_math_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find cmath or math.h.  ml requires at least one or the other.]])
  fi
else
  # TFLOP
  if test "X$ac_cv_header_iomanip" != "Xyes" && test "X$ac_cv_header_iomanip_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find iomanip or iomanip.h.  ml requires at least one or the other.]])
  fi
fi  

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE

# ------------------------------------------------------------------------
# Checks for compiler characteristics - Epetra_ConfigDefs.h uses these, 
# so it makes sense to leave them in - see PORTABILITY in this directory
# ------------------------------------------------------------------------
AC_CXX_NAMESPACES
AC_CXX_HAVE_STL
AC_CXX_BOOL
AC_CXX_MUTABLE
# ---------------------------------------------------------------------
# These are pretty obscure
# ---------------------------------------------------------------------
AC_CXX_NEW_FOR_SCOPING
#AC_CXX_STD_IOS_BASE_FMTFLAGS
AC_CXX_STD_SPRINTF

# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# #np# newp_swahili is an example of a package feature which requires both
# a cpp #define command and one or more source or header files  - ie 
# conditional compilation which requires a change to the makefile.
# If configure is called with "--enable-newp_swahili", it will place 
# a #define HAVE_NEWP_SWAHILI in New_Package_config.h and add source 
# and header files as specified in Makefile.am 
#
# To add an option to your package, you will need to modify this file
# (i.e. configure.ac), and src/Makefile.am.  A bootstrap is also required.
# Look for SWAHILI in both files listed above.
#
# TAC_ARG_ENABLE_OPTION performs several tasks, including:
#   1)  Setting a help string
#   2)  Printing a message to standard out
#   3)  Defining HAVE_NEWP_SWAHILI (For use in New_Package_config.h)
#
# AM_CONDITIONAL makes HAVE_NEWP_SWAHILI available for use in Makefile.am
# ------------------------------------------------------------------------
#
#JJHnp try this option with ML_WITH_EPETRA
#TAC_ARG_ENABLE_OPTION(newp_swahili, [enable swahili support], NEWP_SWAHILI, no)
#AM_CONDITIONAL(HAVE_NEWP_SWAHILI, [test "X$ac_cv_use_newp_swahili" != "Xno"])
# TAC_ARG_WITH_PACKAGE(ml_epetra, [Enable Aztec support], EPETRA, yes)


# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------

#
#  #np# - These options are specific to Epetra and are left only as examples.
#         If uncommented, the following lines would define HAVE_FATAL_MESSAGES
#         in the (package)_config.h file.  Developers should use a 
#         more specific name such as 'new_package-fatal_messages', rather
#         than 'fatal_messages'.
#TAC_ARG_ENABLE_OPTION(fatal_messages, [enable printing of fatal error messages], 
#                   FATAL_MESSAGES, yes)
#TAC_ARG_ENABLE_OPTION(warning_messages, [enable printing of warning error messages], 
#                   WARNING_MESSAGES, no)
#TAC_ARG_ENABLE_OPTION(format_io, [Use package-defined IO formatting for output], 
#                   FORMAT_IO, yes)
#  #np# - These are specific to Epetra and are left only as examples.
#         If uncommented, the following line would define HAVE_ZOLTAN in
#         the (package)_config.h file.  Developers should use a more
#         specific name such as 'ml-zoltan', rather than 'zoltan'.
# TAC_ARG_WITH_PACKAGE(zoltan, [Enable Zoltan interface support], ZOLTAN, no)
# TAC_ARG_WITH_PACKAGE(rtop, [Enable RTOp interface support], RTOP, no)
# TAC_ARG_WITH_PACKAGE(esi, [Enable ESI interface support], ESI, yes)
# TAC_ARG_WITH_PACKAGE(superlu, [Enable SuperLU interface support], SUPERLU, no)

# invoke with "with-xxxxxx", where xxxxxx is the first argument, e.g., ml_aztec

###########################
# these are ON by default
###########################
 TAC_ARG_WITH_PACKAGE(ml_epetra, [Enable Epetra support], ML_EPETRA, yes)
#
 TAC_ARG_WITH_PACKAGE(ml_aztecoo, [Enable Aztecoo support], ML_AZTEC, yes)
#if we're not using aztecoo, then build certain ML examples
#without dependencies on libepetra.a  or libaztecoo.a
#(See ml/examples/Makefile.am.)
AM_CONDITIONAL(ML_STANDALONE_MODE, test "X$ac_cv_use_ml_aztecoo" != "Xyes")
# TAC_ARG_WITH_PACKAGE(ml_mpi, [Enable mpi support], ML_MPI, yes)
 TAC_ARG_WITH_PACKAGE(ml_external_mpi_functions, [Use mpi-provided global communication functions (not ml's internal implementation)], ML_EXTERNAL_MPI_FUNCTIONS, yes)
# TAC_ARG_WITH_PACKAGE(ml_vendor_blas, [Use system blas library (not ml's internal reference blas)], ML_VENDOR_BLAS, yes)
# TAC_ARG_WITH_PACKAGE(ml_vendor_lapack, [Use system lapack library (not ml's internal reference lapack)], ML_VENDOR_LAPACK, yes)
 TAC_ARG_ENABLE_OPTION(ml_enrich, [Enable latest enriched Maxwell prolongator], ML_ENRICH, no)
#
 TAC_ARG_ENABLE_OPTION(ml_newtpe, [Enable latest Maxwell prolongator fixes], ML_NEW_T_PE, yes)

###########################
# these are OFF by default
###########################
 TAC_ARG_WITH_PACKAGE(ml_superlu, [Enable SuperLU support.  Note 1: SuperLU is a third party library.  As such, configure assumes that it has already been built.  Note 2: You must also specify the location and name of the SuperLU headers and libraries like so: configure --with-libs="<other libraries> -lsuperlu" --with-ldflags="-L<SuperLU-library-path>" --with-incdirs="<SuperLU-include-path>".], ML_SUPERLU, no)
#
 TAC_ARG_WITH_PACKAGE(ml_superlu_dist, [Enable SuperLU_DIST support.  Note 1: SuperLU_DIST is a third party library.  As such, configure assumes that it has already been built.  Note 2: You must also specify the location and name of the SuperLU_DIST headers and libraries like so: configure --with-libs="<other libraries> -lsuperlu_dist" --with-ldflags="-L<SuperLU_DIST-library-path>" --with-incdirs="<SuperLU_DIST-include-path>".], ML_SUPERLU_DIST, no)
 
# Note 2: You must also specify the location of the SuperLU_DIST headers and libraries via the configure options --with-ldflags=<SuperLU_DIST-library-path/library-name> and --with-incdirs=SuperLU_DIST-include-directory.], ML_SUPERLU_DIST, no)
#
 TAC_ARG_WITH_PACKAGE(ml_aztec2_1, [Enable Aztec 2.1 support], ML_AZTEC2_1, no)
#
 TAC_ARG_ENABLE_OPTION(ml_complex_maxwell, [Enable complex maxwell solver capability], ML_COMPLEX_MAXWELL, no)
 #
#
 TAC_ARG_ENABLE_OPTION(ml_timing, [Enable timing routines], ML_TIMING, no)
#
 TAC_ARG_ENABLE_OPTION(ml_flops, [Enable flop counting.], ML_FLOPS, no)
#
 TAC_ARG_ENABLE_OPTION(ml_benchmarking, [For benchmarking suite.], ML_BENCHMARKING, no)
AM_CONDITIONAL(ML_BENCHMARKING, test "X$ac_cv_use_ml_benchmarking" != "Xno")

#TAC_ARG_ENABLE_FEATURE_SUB( ml, slud, [Enable ML interface support for SuperLU], ML_SLU, no)

#if test "X$ac_cv_header_cstdlib" != "Xyes" && test "X$ac_cv_header_stdlib_h" != "Xyes"; then
#  AC_MSG_ERROR([[Cannot find cstdlib or stdlib.h.  ml requires at least one or the other.]])
#fi

TAC_ARG_WITH_PACKAGE(ml_metis, [Enable METIS support.  METIS can be used to create aggregates using graph-decomposition algorithsm. Note 1: METIS is a third party library.  As such, configure assumes that it has already been built.  Note 2: You must also specify the location of headers and libraries like, for instance: configure --with-libs="<other libraries> -lmetis" --with-ldflags="-L<METIS-library-path>" and --with-incdirs="<METIS-include-path>".], ML_METIS, no)
#
TAC_ARG_WITH_PACKAGE(ml_parmetis2x, [Enable ParMETIS support, version 2.x. See notes for --with-ml_metis.], ML_PARMETIS_2x, no) 
#
TAC_ARG_WITH_PACKAGE(ml_parmetis3x, [Enable ParMETIS support, version 3.x. See notes for --with-ml_metis.], ML_PARMETIS_3x, no) 
#
TAC_ARG_WITH_PACKAGE(ml_amesos, [Enable AMESOS support. Amesos, a Trilinos package, can be used to interface ML with various direct solvers, such as SuperLU, MUMPS, UMFPACK, KLU, and other. Amesos requires --with-ml_epetra. See notes for --enable-ml_metis.], ML_AMESOS, no)
AM_CONDITIONAL(ML_USING_AMESOS, test "X$ac_cv_use_ml_amesos" = "Xyes")
#
TAC_ARG_WITH_PACKAGE(ml_ifpack, [Enable IFPACK support. Ifpack, a Trilinos package, can provide incomplete factorization smoothers to ML. Ifpack requires --with-ml_epetra. See notes for --enable-ml_metis.], ML_IFPACK, no)
AM_CONDITIONAL(ML_USING_IFPACK, test "X$ac_cv_use_ml_ifpack" = "Xyes")
#
TAC_ARG_WITH_PACKAGE(ml_teuchos, [Enable TEUCHOS support. Teuchos, a Trilinos package, is required by the class Epetra_ML_Preconditioner. Teuchos requires --with-ml_epetra. See notes for --enable-ml_metis.], ML_TEUCHOS, no)
AM_CONDITIONAL(ML_USING_TEUCHOS, test "X$ac_cv_use_ml_teuchos" = "Xyes")
#
TAC_ARG_WITH_PACKAGE(ml_anasazi, [Enable ANASAZI support. Anasazi, a Trilinos package, can be used for eigenvalue estimates. Anasazi requires --with-ml_epetra and --with-ml_teuchos.], ML_ANASAZI, no)
AM_CONDITIONAL(ML_USING_ANASAZI, test "X$ac_cv_use_ml_anasazi" = "Xyes")
#
TAC_ARG_WITH_PACKAGE(ml_triutils, [Enable TRIUTILS support. Triutils, a Trilinos package, is used if some ML examples to construct the linear system. Triutils requires --with-ml_epetra.], ML_TRIUTILS, no)
AM_CONDITIONAL(ML_USING_TRIUTILS, test "X$ac_cv_use_ml_triutils" = "Xyes")
#
AC_ARG_WITH(ml-64bit-integer,
[AC_HELP_STRING([--with-ml-64bit-integer],[Support for 64-bit integers; required if ML used
with MPSalsa.  Examples: --with-ml-64bit-integer="long long" for LINUX, TFLOP, IRIXN32"; --with-ml-64bit-integer="long" for IRIX64, ALPHA (Cplant)])],
[ML_BIG_INT_SUPPORT=${withval}],
[ML_BIG_INT_SUPPORT="int"]
)
AC_DEFINE_UNQUOTED(ML_BIG_INT, $ML_BIG_INT_SUPPORT,
   [Support for 64-bit integers])

if test "X${ML_BIG_INT_SUPPORT}" = "Xyes" && test "X${ML_BIG_INT}" = "X" || test "X${ML_BIG_INT_SUPPORT}" = "X" && test "X${ML_BIG_INT}" = "X"; then
  echo "ERROR: An argument to --with-ml-64bit-integer must be specified!"
  echo "  For example:"
  echo "    --with-ml-64bit-integer=\"long long\" for LINUX, TFLOP, IRIXN32"
  echo "    --with-ml-64bit-integer=\"long\"      for IRIX64, ALPHA (Cplant)"
  exit 1
fi

# ------------------------------------------------------------------------
# Checks if tests and examples should be built
# ------------------------------------------------------------------------

#TAC_ARG_ENABLE_FEATURE(tests, [Build tests for all Trilinos packages (not all packages are sensitive to this option)], TESTS, yes)
#AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_tests" != "Xno")

TAC_ARG_ENABLE_FEATURE(examples, [Build examples for ML package], EXAMPLES, yes)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_examples" != "Xno")

# ------------------------------------------------------------------------
# Specify other directories 
# ------------------------------------------------------------------------

# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS

# Yet another opportunity to remove code if you aren't
# using Fortran
# Define F77_FUNC that will be used to link with Fortran subroutines, e.g., BLAS/LAPACK - trash WORKGXX 
AC_F77_WRAPPERS

# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

# #np# - 
# To add any libraries that you want to check for, replace pow below 
# with the name of an entry point in that library.
# checks for libraries now
dnl Replace `main' with a function in -lm:
#  #np# These lines can be removed if the math lib is not needed.
AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))

TAC_ARG_WITH_BLASLIB
TAC_ARG_WITH_LAPACKLIB

# find or choose to supply valid blas and lapack libraries
  ACX_BLAS([], [AC_MSG_ERROR([[Cannot find blas library, specify a path using --with-ldflags=-I<DIR> (ex. --with-ldflags=-I/usr/path/lib) or a specific library using --with-blas=DIR/LIB (for example --with-blas=/usr/path/lib/libcxml.a)]])])
  ACX_LAPACK([], [AC_MSG_ERROR([[Cannot find lapack library, specify a path using --with-ldflags=-I<DIR> (ex. --with-dlflags=-I/usr/path/lib) or a specific library using --with-lapack=DIR/LIB (ex.  --with-lapack=/usr/path/lib/libcxml.a)]])])

# ------------------------------------------------------------------------
# Checks for compiler characteristics
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for linker characteristics
# ------------------------------------------------------------------------

# Determine libraries needed for linking with Fortran
AC_F77_LIBRARY_LDFLAGS


#  ----------------------------------------------------
#  Needs this because Dec handles templates differently
#  ----------------------------------------------------
#  Only uncomment the following lines if the package
#  creates a subdirectory called cxx_repository for
#  template object files on DEC.  (HINT: try to build on 
#  DEC/CPLANT and see if template-related errors occur)

# -- DEC CXX --
#case $target in
#  alpha*)
#    if test $CXX = cxx; then
#      TRILINOS_TEMPLATE_OBJS="cxx_repository/*.o"
#    fi
#  ;;
#esac

#AC_MSG_CHECKING([for template-specific object files to be linked])
#if test "X$TRILINOS_TEMPLATE_OBJS" = "X"; then
# AC_MSG_RESULT(none)
#else
# AC_MSG_RESULT($TRILINOS_TEMPLATE_OBJS)
#fi

#AC_SUBST(TRILINOS_TEMPLATE_OBJS)

# ------------------------------------------------------------------------
# Checks for library functions
# ------------------------------------------------------------------------

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
#  Most packages should have an example and a full suite 
#  of tests.  New package comes only with src and example.
#  You will need to change AC_CONFIG_FILES below and Makefile.am
#  to add a new directory.
AC_CONFIG_FILES([
		Makefile 
		src/Makefile
		examples/Makefile
		])

AC_OUTPUT()

# Bye World!
echo "----------------------------------------"
echo "Finished Running ml Configure Script"
echo "----------------------------------------"
