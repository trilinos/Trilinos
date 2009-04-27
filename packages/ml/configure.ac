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
AC_INIT(ml, 4.0, jhu@sandia.gov)

# Hello World!
#
echo "----------------------------------------"
echo "Running ML Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.  Any file in that directory will do,
# preferably one that is unlikely to be removed or renamed.

AC_CONFIG_SRCDIR([src/Main/ml_struct.c])

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

# Use automake

AM_INIT_AUTOMAKE(1.10 no-define tar-ustar)

# Specify required version of autoconf.

AC_PREREQ(2.61)

#  can eliminate compiler checks below if your package does not use the
#  language corresponding to the check
# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI

AC_PROG_CXX(CC g++ c++ cxx)

AC_ARG_ENABLE( ml_multiple_rhs,
             AC_HELP_STRING([--enable-ml_multiple_rhs],[Enables blocked code
                for multiple right hand sides.  Example:
                --enable-ml_multiple_rhs=4 for a blocking factor of 4.
                This option requires that Epetra be enable.
                NOTE:  With the icc compiler, you may need to add
                CXXFLAGS="-Kc++" to the configure line, in order for it to
                interpret *.c files as c++ files.  Similarly, with Compaq
                compilers, you may need to add CXXFLAGS="-x cxx" to the
                configure line.  One word of warning:  be sure that whatever
                compiler flags you want used are given in CXXFLAGS, as CFLAGS
                will be overwritten by CXXFLAGS. (default is no)]),
             [ac_cv_use_ml_multiple_rhs=${enableval}],
             [ac_cv_use_ml_multiple_rhs=no])

if test "X$ac_cv_use_epetra" == "Xno" && test "X$ac_cv_use_ml_multiple_rhs" != "Xno"; then
 AC_MSG_ERROR([[The use of --enable-ml_multiple_rhs in ML requires Epetra.]])
# ac_cv_use_ml_multiple_rhs=no;
fi

if test "X$ac_cv_use_ml_multiple_rhs" != "Xno"; then
CC=${CXX}
# if multiple rhs are going to be used, create a preprocessor macro
#  that has the value given on the configure line
AC_DEFINE_UNQUOTED(ML_MULTIPLE_RHS_BLOCK_FACTOR, $ac_cv_use_ml_multiple_rhs,
   [Support for multiple right hand sides.])
fi

AC_PROG_CC(cc gcc)
AC_PROG_F77(f77 g77 gfortran f90 xlf90 f95)
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Get the absolute path to the build directory.
#
# We need this to be able to figure out paths to the libraries for the
# other Trilinos packages.
# ------------------------------------------------------------------------

AC_MSG_CHECKING([finding build directory])
ML_BUILD_DIRECTORY=`pwd`
AC_MSG_RESULT($ML_BUILD_DIRECTORY)

# ------------------------------------------------------------------------
# Perform substitutions in output files
# ------------------------------------------------------------------------
                                                                                
##AC_SUBST(ac_aux_dir)
AC_SUBST(ML_BUILD_DIRECTORY)
                                                                                
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
# The next line is to override the CFLAGS line.  This is necessary for
# --with-ml_multiple_rhs to work properly on stratus.  Bsically, the problem is
# that the cxx compiler needs an additional option "-x cxx" to force it to
# compile all sources as if they are c++ sources.
# 
# configure uses cc and CFLAGS during its internal builds.  On stratus the cc
# compiler doesn't understand "-x cxx", so this option can't be set in CFLAGS
# on the configure line.
# Just setting "-x cxx" via CXXFLAGS on the configure line doesn't work because
# cxx in the end will be using the options in CFLAGS.
#
# The solution is to set CXXFLAGS = "-x cxx" and then override the value of
# CFLAGS with CXXFLAGS after configure compile testing is done. 
#
# This is a very fragile option, one I view as a temporary patch until some of
# the multiple right hand side code is rewritten.
#
if test "X$ac_cv_use_ml_multiple_rhs" != "Xno"; then
CFLAGS=${CXXFLAGS}
fi
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

# ------------------------------------------------------------------------
# Checks for Makefile.export related systems
# ------------------------------------------------------------------------
TAC_ARG_ENABLE_EXPORT_MAKEFILES(yes)

# Customize the checks for headers needed in package being added.
# NOTE: This section still needs to be updated for ml, it is still 
# code from Epetra
# ------------------------------------------------------------------------
# Checks for header files - used in Epetra_ConfigDefs.h 
# ------------------------------------------------------------------------

AC_CHECK_HEADERS([cstdlib stdlib.h], [break],)
AC_CHECK_HEADERS([cstdio stdio.h], [break],)
AC_CHECK_HEADERS([cassert assert.h], [break],)
AC_CHECK_HEADERS([cmath math.h], [break],)
AC_CHECK_HEADERS([iomanip iomanip.h], [break],)
AC_CHECK_HEADERS([iostream iostream.h], [break],)
AC_CHECK_HEADERS([string string.h], [break],)
AC_CHECK_HEADERS([malloc.h],,)

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

if test "X$ac_cv_header_malloc_h" != "Xyes"; then
 AC_FUNC_MALLOC
 #if test "X$ac_cv_func_malloc_works" != "Xyes"; then
 #  AC_MSG_ERROR([malloc function not working, needed by TriUtils.])
 #fi
fi

# Checks for typedefs, structures, and compiler characteristics.
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

AC_MSG_CHECKING([whether $CC recognizes __PRETTY_FUNCTION__])
AC_TRY_COMPILE(, [const char *s = __PRETTY_FUNCTION__; ],
   AC_DEFINE(HAVE_CFUNC,, 
   [Define if the C complier supports __PRETTY_FUNCTION__]) ac_cv_c_var_func=yes,
   ac_cv_c_var_func=no)
AM_CONDITIONAL(ML_CFUNC, test "X$ac_cv_c_var_func" != "Xno")

if test "X$ac_cv_c_var_func" != "Xno"; then
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
fi



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

#####################################
# Following options are ON by default
#####################################

TAC_ARG_ENABLE_OPTION(ml_newtpe, [Enable latest Maxwell prolongator fixes], ML_NEW_T_PE, yes)

#######################################
# Following options are OFF by default
#######################################
 TAC_ARG_ENABLE_OPTION(ml_enrich, [Enable latest enriched Maxwell prolongator], ML_ENRICH, no)
#
 TAC_ARG_WITH_PACKAGE(ml_aztec2_1, [Enable Aztec 2.1 support], ML_AZTEC2_1, no)
#
 TAC_ARG_ENABLE_OPTION(ml_complex_maxwell, [Enable complex maxwell solver capability], ML_COMPLEX_MAXWELL, no)
 #
#
 TAC_ARG_ENABLE_OPTION(ml_timing, [Enable timing routines], ML_TIMING, no)
#
 TAC_ARG_ENABLE_OPTION(ml_memory_checking, [Enable ml memory checking.],
ML_MEMORY_CHECK, no)
#
 TAC_ARG_ENABLE_OPTION(ml_flops, [Enable flop counting.], ML_FLOPS, no)
#
 TAC_ARG_ENABLE_OPTION(ml_benchmarking, [For benchmarking suite.], ML_BENCHMARKING, no)
AM_CONDITIONAL(ML_BENCHMARKING, test "X$ac_cv_use_ml_benchmarking" != "Xno")
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
#
############################################
#### interfaces to 3rd party packages
############################################
#
TAC_ARG_WITH_PACKAGE(ml_superlu, [(Deprecated in favor of --enable-ml-amesos.
Cannot be used at the same time as --enable-ml-amesos.) Enable serial SuperLU (version 1.0) support.  Note 1: SuperLU is a third party library.  As such, configure assumes that it has already been built.  Note 2: You must also specify the location and name of the SuperLU headers and libraries like so: configure --with-libs="<other libraries> -lsuperlu" --with-ldflags="-L<SuperLU-library-path>" --with-incdirs="<SuperLU-include-path>".], ML_SUPERLU, no)
#
TAC_ARG_WITH_PACKAGE(ml_superlu2, [(Deprecated in favor of --enable-ml-amesos.  Cannot be used at the same time as --enable-ml-amesos.) Enable serial SuperLU (version 2.0) support.  Note 1: SuperLU is a third party library.  As such, configure assumes that it has already been built.  Note 2: You must also specify the location and name of the SuperLU headers and libraries like so: configure --with-libs="<other libraries> -lsuperlu" --with-ldflags="-L<SuperLU-library-path>" --with-incdirs="<SuperLU-include-path>".], ML_SUPERLU2_0, no)
#
TAC_ARG_WITH_PACKAGE(ml_superlu_dist, [(Deprecated in favor of --enable-ml-amesos.  Cannot be used at the same time as --enable-ml-amesos.) Enable SuperLU_DIST support.  Note 1: SuperLU_DIST is a third party library.  As such, configure assumes that it has already been built.  Note 2: You must also specify the location and name of the SuperLU_DIST headers and libraries like so: configure --with-libs="<other libraries> -lsuperlu_dist" --with-ldflags="-L<SuperLU_DIST-library-path>" --with-incdirs="<SuperLU_DIST-include-path>".], ML_SUPERLU_DIST, no)
#
TAC_ARG_WITH_PACKAGE(ml_metis, [Enable METIS support.  METIS can be used to create aggregates using graph-decomposition algorithsm. Note 1: METIS is a third party library.  As such, configure assumes that it has already been built.  Note 2: You must also specify the location of headers and libraries like, for instance: configure --with-libs="<other libraries> -lmetis" --with-ldflags="-L<METIS-library-path>" and --with-incdirs="<METIS-include-path>".], ML_METIS, no)
#
TAC_ARG_WITH_PACKAGE(ml_parmetis, [Enable ParMETIS support. See notes for --with-ml_metis.], ML_PARMETIS, no) 
#
TAC_ARG_WITH_PACKAGE(ml_arpack, [Enable ARPACK support. See notes for --with-ml_metis.], ML_ARPACK, no) 
#
TAC_ARG_WITH_PACKAGE(ml_parpack, [Enable PARPACK support. See notes for --with-ml_metis.], ML_PARPACK, no) 
#
TAC_ARG_WITH_PACKAGE(ml_parasails, [Enable Parasails support. See notes for --with-ml_metis.], ML_PARASAILS, no) 
#
TAC_ARG_WITH_PACKAGE(ml_nox, [Enable ML_NOX support. --with-ml_nox requires that configure and make with --enable-nox --enable-prerelease --enable-nox-epetra --enable-epetraext have been run before, otherwise it won't find necessary headers and won't link. ], ML_NOX, no) 
AM_CONDITIONAL(ML_USING_NOX, test "X$ac_cv_use_ml_nox" = "Xyes")
#
TAC_ARG_WITH_PACKAGE(ml_zoltan, [Enable Zoltan support. This option requires that MPI be enabled.  See notes for --with-ml_metis.], ML_ZOLTAN, no) 
AM_CONDITIONAL(ML_USING_ZOLTAN, test "X$ac_cv_use_ml_zoltan" = "Xyes")
#
TAC_ARG_WITH_PACKAGE(ml_zoltan3, [Enable Zoltan 3.0 support. This option requires that MPI and --enable-zoltan.], ML_ZOLTAN_THREE, no)
#
if test "X${HAVE_PKG_MPI}" != "Xyes"; then
   #issue an error message and abort if the user has configured with
   #the zoltan interface enabled but without mpi enabled
   if test "X$ac_cv_use_ml_zoltan" = "Xyes" || test "X$ac_cv_use_ml_zoltan3" = "Xyes"; then
     AC_MSG_ERROR([[Zoltan support in ML requires enabling MPI.]])
   fi
fi
#
TAC_ARG_WITH_PACKAGE(petsc, [Enable PETSc support.], PETSC, no)

#
# ------------------------------------------------------------------------
# Interfaces to other Trilinos packages
# ------------------------------------------------------------------------
#
TAC_ARG_ENABLE_CAN_USE_PACKAGE(ml, epetra, ML_EPETRA,  ML_USING_EPETRA, yes, [Enable Epetra support.], [Enable Epetra support. Can be overridden with --disable-ml-epetra.])
#
TAC_ARG_ENABLE_CAN_USE_PACKAGE(ml, amesos, ML_AMESOS, ML_USING_AMESOS, yes, [Enable AMESOS support. Amesos, a Trilinos package, can be used to interface ML with various direct solvers, such as SuperLU, MUMPS, UMFPACK, KLU, and other.  Amesos requires --enable-epetra. See notes for --enable-ml_metis.], [Enable AMESOS support. Amesos, a Trilinos package, can be used to interface ML with various direct solvers, such as SuperLU, SuperLUdist, MUMPS, UMFPACK, KLU, and other. Amesos requires --enable-epetra. See notes for --enable-ml_metis. Can be overridden with --disable-ml-amesos.])
#
# JJH 7/27/2006 commenting this next section out because it causes problems for
# ALEGRA
#
#if test "X${ac_cv_use_amesos}" = "Xyes"; then
#   #Issue an error message and abort if the user has configured with
#   #both the Amesos interface and ML's internal SuperLU interface enabled.
#   #Amesos wants SuperLU version 3.0, whereas ML uses SuperLU versions 1.0/2.0.
#   if test "X$ac_cv_use_ml_superlu" = "Xyes" || test "X${ac_cv_use_ml_superlu2}" = "Xyes" || test "X${ac_cv_use_ml_superlu_dist}" = "Xyes"; then
#     AC_MSG_ERROR([[You cannot use Amesos and ML's internal interface to SuperLU/SuperLU-dist at the same time.]])
#   fi
#fi
#
TAC_ARG_ENABLE_CAN_USE_PACKAGE(ml, ifpack, ML_IFPACK, ML_USING_IFPACK, yes, [Enable IFPACK support. Ifpack, a Trilinos package, can provide incomplete factorization smoothers to ML. Ifpack requires --enable-epetra. See notes for --enable-ml_metis.], [Enable IFPACK support. Ifpack, a Trilinos package, can provide incomplete factorization smoothers to ML. Ifpack requires --enable-epetra. See notes for --enable-ml_metis. Can be overridden with --disable-ml-ifpack.])
#
TAC_ARG_ENABLE_CAN_USE_PACKAGE(ml, teuchos, ML_TEUCHOS, ML_USING_TEUCHOS, yes, [Enable TEUCHOS support. Teuchos, a Trilinos package, is required by the class Epetra_ML_Preconditioner. Teuchos requires --enable-epetra. See notes for --enable-ml_metis.], [Enable TEUCHOS support. Teuchos, a Trilinos package, is required by the class Epetra_ML_Preconditioner. Teuchos requires --enable-epetra. See notes for --enable-ml_metis. Can be overridden with --disable-ml-teuchos.])
#
#TAC_ARG_ENABLE_CAN_USE_PACKAGE(ml, anasazi, ML_ANASAZI, ML_USING_ANASAZI, yes, [Enable ANASAZI support. Anasazi, a Trilinos package, can be used for eigenvalue estimates. Anasazi requires --enable-epetra and --enable-teuchos.], [Enable ANASAZI support. Anasazi, a Trilinos package, can be used for eigenvalue estimates. Anasazi requires --enable-epetra and --enable-teuchos. Can be overridden with --disable-ml-anasazi.])
#
TAC_ARG_ENABLE_CAN_USE_PACKAGE(ml, aztecoo, ML_AZTECOO, ML_USING_AZTECOO, yes, [Enable AztecOO support. AztecOO, a Trilinos package, can be used in adaptive preconditioning.  Requires --enable-epetra --enable-teuchos. See notes for --enable-ml_metis.], [Enable AztecOO support. AztecOO, a Trilinos package, can be used in adaptive preconditioning.  Requires --enable-epetra --enable-teuchos. See notes for --enable-ml_metis. Can be overridden with --disable-ml-aztecoo.])
#
TAC_ARG_ENABLE_CAN_USE_PACKAGE(ml, epetraext, ML_EPETRAEXT, ML_USING_EPETRAEXT, yes, [Enable EpetraExt support. EpetraExt requires --enable-epetra. See notes for --enable-ml_metis.], [Enable EpetraExt support. EpetraExt requires --enable-epetra. See notes for --enable-ml_metis. Can be overridden with --disable-ml-epetraext.])
#
# if user wants PETSc, EpetraExt must be enabled.
if test "X$ac_cv_use_epetraext" != "Xyes"; then
   if test "X$ac_cv_use_petsc" = "Xyes"; then
     AC_MSG_ERROR([[--with-petsc requires --enable-epetraext.]])
   fi
fi
#
TAC_ARG_ENABLE_CAN_USE_PACKAGE(ml, isorropia, ML_ISORROPIA, ML_USING_ISORROPIA, yes, [Enable Isorropia support. Isorropia requires --enable-isorropia. See notes for --enable-ml_metis.], [Enable Isorropia support. Isorropia requires --enable-isorropia. See notes for --enable-ml_metis. Can be overridden with --disable-ml-isorropia.])
#
TAC_ARG_ENABLE_CAN_USE_PACKAGE(ml, galeri, ML_GALERI, ML_USING_GALERI, yes, [Enable Galeri support. Galeri, a Trilinos package, is used in some ML examples to construct the linear system and can optionally be used by the MLAPI.], [Enable Galeri support.  Galeri, a Trilinos package, is used in some ML examples to construct the linear system and can optionally be used by the MLAPI.  Can be overridden with --disable-ml-galeri.])

#if we're not using epetra, then build certain ML examples
#without dependencies on libepetra.a (and use C linker for C files,
# otherwise use C++ linker for C examples as well)
#(See ml/examples/Makefile.am.)
AM_CONDITIONAL(ML_STANDALONE_MODE, test "X$ac_cv_use_epetra" != "Xyes")

# ------------------------------------------------------------------------
# Checks if tests and examples should be built
# ------------------------------------------------------------------------

# Embedded version of tac_arg_enable_feature_sub_check.m4 to support
# deprecated --enable-ml_examples and also --enable-ml-examples

AC_ARG_ENABLE(examples, AC_HELP_STRING([--enable-examples], [Make examples for all Trilinos packages buildable with 'make examples' (default is yes)]), ac_cv_use_examples=$enableval, ac_cv_use_examples=yes)

AC_ARG_ENABLE(ml_examples, AC_HELP_STRING([--enable-ml_examples], [Make ML examples buildable with 'make examples' (deprecated, use --enable-ml-examples)]), ac_cvs_use_ml_examples=$enableval, ac_cv_use_ml_examples=${ac_cv_use_examples})

AC_ARG_ENABLE([ml-examples],
AC_HELP_STRING([--enable-ml-examples],[Make ML examples buildable with 'make examples' (default is yes if --disable-examples is not specified)]),
ac_cv_use_ml_examples=$enableval, )

AC_MSG_CHECKING(whether to use ml-examples)

if test "X$ac_cv_use_ml_examples" != "Xno"; then
  AC_MSG_RESULT(yes)
  AC_DEFINE([HAVE_ML_EXAMPLES],,[Define if want to build ml-examples])
else
  AC_MSG_RESULT(no)
fi

AM_CONDITIONAL(ML_BUILD_EXAMPLES, test "X$ac_cv_use_ml_examples" != "Xno")

# End of embedded tac_arg_enable_feature_sub_check.m4

TAC_ARG_ENABLE_FEATURE(tests, [Make tests for all Trilinos packages buildable with 'make tests'], TESTS, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( ml, tests, [Make ML tests buildable with 'make tests'], ML_TESTS)
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_ml_tests" != "Xno")

#We now build tests and examples through separate make targets, rather than
#during "make".  We still need to conditionally include the test and example
#in SUBDIRS, even though SUB_TEST and SUB_EXAMPLE will never be
#defined, so that the tests and examples are included in the distribution
#tarball.
AM_CONDITIONAL(SUB_TEST, test "X$ac_cv_use_sub_test" = "Xyes")
AM_CONDITIONAL(SUB_EXAMPLE, test "X$ac_cv_use_sub_example" = "Xyes")

TAC_ARG_ENABLE_FEATURE(libcheck, [Check for some third-party libraries including BLAS and LAPACK.  (Cannot be disabled unless tests and examples are also disabled.)], LIBCHECK, yes)

# ------------------------------------------------------------------------
# Specify other directories 
# ------------------------------------------------------------------------

# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS


# Yet another opportunity to remove code if you aren't
# using Fortran
# Define F77_FUNC that will be used to link with Fortran subroutines.
AC_F77_WRAPPERS

# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

# If tests, examples and libcheck are disabled, we don't have to check
# for these libraries.
if test "X$ac_cv_use_ml_tests" != "Xno" || test "X$ac_cv_use_ml_examples" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then

  # #np# - 
  # To add any libraries that you want to check for, replace pow below 
  # with the name of an entry point in that library.
  # checks for libraries now
  dnl Replace `main' with a function in -lm:
  #  #np# These lines can be removed if the math lib is not needed.
  AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
  AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))

  # find or choose to supply valid blas and lapack libraries
  TAC_BLAS([], [AC_MSG_ERROR([[Cannot find blas library, specify a path using --with-ldflags=-I<DIR> (ex. --with-ldflags=-I/usr/path/lib) or a specific library using --with-blas=DIR/LIB (for example --with-blas=/usr/path/lib/libcxml.a)]])])
  TAC_LAPACK([], [AC_MSG_ERROR([[Cannot find lapack library, specify a path using --with-ldflags=-I<DIR> (ex. --with-ldflags=-I/usr/path/lib) or a specific library using --with-lapack=DIR/LIB (ex.  --with-lapack=/usr/path/lib/libcxml.a)]])])

else

  # If we don't check for BLAS above, we will use F77_FUNC as the name mangling
  # for the BLAS and LAPACK.  This is what we used before we allowed the BLAS
  # and LAPACK to use different name mangling than the compiler.
  AC_F77_FUNC(dgemm)
  if test "X$dgemm" = "Xdgemm_"; then
    AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [name ## _], [Define the Fortran name mangling to be used for the BLAS])
  else
    if test "X$dgemm" = "Xdgemm"; then
      AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [name], [Define the Fortran name mangling to be used for the BLAS])
    else
      if test "X$dgemm" = "XDGEMM"; then
        AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [NAME], [Define the Fortran name mangling to be used for the BLAS])
      else
        if test "X$dgemm" = "XDGEMM_"; then
           AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [NAME ## _], [Define the Fortran name mangling to be used for the BLAS])
        else
          AC_MSG_ERROR([[Cannot determine BLAS and LAPACK name mangling using the value of F77_FUNC.  Please try configuring Trilinos without the --disable-libcheck option.]])
        fi
      fi
    fi
  fi

fi
# end of the list of libraries that don't need to be checked for if
# tests and examples are disabled.


# ------------------------------------------------------------------------
# Enable optional matlab wrappers - this needs to be after the other directories
# or else we'll never find the matlab libraries / includes
# ------------------------------------------------------------------------
TAC_ARG_ENABLE_FEATURE_SUB( ml, matlab, [Build ML Matlab wrappers (requires --with-matlab, which is NOT set by default)], ML_MATLAB, no)
AM_CONDITIONAL(BUILD_ML_MATLAB, test "X$ac_cv_use_ml_matlab" != "Xno")
TAC_ENABLE_MATLAB

if test "X${HAVE_PKG_MPI}" != "Xno"; then
   #issue an error message and abort if the user has configured with
   # Matlab and MPI
   if test "X$ac_cv_use_ml_matlab" = "Xyes"; then
     AC_MSG_ERROR([[MLMEX requires disabling MPI.]])
   fi
fi



# ------------------------------------------------------------------------
# Checks for compiler characteristics
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for linker characteristics
# ------------------------------------------------------------------------

AC_MSG_CHECKING(for mallinfo() function)
AC_TRY_LINK_FUNC(mallinfo,,mlmallinfook="yes",mlmallinfook="no")

if test "X$mlmallinfook" != "Xyes"; then
  AC_DEFINE(ML_MALLINFO,[],[mallinfo with ML])
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
fi

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
if test "X$TRILINOS_TEMPLATE_OBJS" = "X"; then
 case $target in
  alpha*)
    if test $CXX = cxx; then
      TRILINOS_TEMPLATE_OBJS="cxx_repository/*.o"
    fi
  ;;
 esac
fi

AC_MSG_CHECKING([for template-specific object files to be linked])
if test "X$TRILINOS_TEMPLATE_OBJS" = "X"; then
 AC_MSG_RESULT(none)
else
 AC_MSG_RESULT($TRILINOS_TEMPLATE_OBJS)
fi

AC_SUBST(TRILINOS_TEMPLATE_OBJS)

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
#AC_CONFIG_FILES([
#		Makefile 
#		Makefile.export
#		src/Makefile
#		examples/Makefile
#                test/Makefile
#                test/MultiLevelPreconditioner_NonSym/Makefile
#		test/MultiLevelPreconditioner_Sym/Makefile
#		test/MultiLevelPreconditioner_Maxwell/Makefile
#		test/ML_Operator2Epetra_RowMatrix/Makefile
#		test/Scalability/Makefile
#		test/VerySimple/Makefile
#		test/AztecSimple/Makefile
#		test/2d_Poisson/Makefile
#		test/3d_Poisson/Makefile
#		test/Zoltan/Makefile
#		test/IFPACKSmoothers/Makefile
#                test/scripts/Makefile
#                test/scripts/daily/Makefile
#                test/scripts/daily/mpi/Makefile
#                test/scripts/daily/serial/Makefile
#		])
AC_CONFIG_FILES([
        Makefile 
        Makefile.export.ml
        src/Makefile
        examples/Makefile
        examples/Advanced/Makefile
        examples/BasicExamples/Makefile
        examples/Elasticity/Makefile
        examples/MatrixFormats/Makefile
        examples/MatrixReading/Makefile
        examples/MatrixFree/Makefile
        examples/StandAlone/Makefile
        examples/TwoLevelDD/Makefile
        examples/Maxwell/Makefile
        examples/Visualization/Makefile
        examples/MultipleRHS/Makefile
        examples/MLAPI/Makefile
        examples/XML/Makefile
        examples/Other/Makefile
        test/Makefile
        test/MultiLevelPreconditioner_NonSym/Makefile
        test/MultiLevelPreconditioner_Sym/Makefile
        test/MultiLevelPreconditioner_Maxwell/Makefile
        test/MultiLevelPreconditioner_Aux/Makefile
        test/MultiLevelPreconditioner_Restart/Makefile
        test/MultiLevelPreconditioner_Reduce/Makefile
        test/ML_Operator2Epetra_RowMatrix/Makefile
        test/AztecSimple/Makefile
        test/2d_Poisson/Makefile
        test/3d_Poisson/Makefile
        test/IFPACKSmoothers/Makefile
        test/MLAPICompare/Makefile
        test/MLAPILinearCombination/Makefile
        test/MLAPISerialMatrix/Makefile
        test/MLAPIDistributedMatrix/Makefile
        test/Communicator/Makefile
        test/SelfSmoother/Makefile
        test/METISCoarsening/Makefile
        test/MatrixFree/Makefile
        test/Miscellaneous/Makefile
        test/RefMaxwell/Makefile
        test/Scalability/Makefile
        test/Zoltan/Makefile
        test/Zoltan_hypergraph/Makefile
        test/ValidateParams/Makefile
        test/BlockCheby/Makefile
        matlab/Makefile
                ])

AC_OUTPUT()

# Bye World!
echo "----------------------------------------"
echo "Finished Running ml Configure Script"
echo "----------------------------------------"
