# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
# 
#                    Teuchos: Common Tools Package
#                 Copyright (2004) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# Initialization 
# ------------------------------------------------------------------------



# This must be the first line in configure.ac.
# Optional 3rd argument is email address for bugs.

#  #np# - package name, version number, and e-mail address below
AC_INIT(teuchos, 1.6d, hkthorn@sandia.gov)

# Hello World!
#  #np# - package name below
echo "----------------------------------------"
echo "Running Teuchos Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.  Any file in that directory will do,
# preferably one that is unlikely to be removed or renamed.

# #np# - file name below
AC_CONFIG_SRCDIR([src/Teuchos_LAPACK.hpp])

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

#
# #np# This a a subdirectory of Trilinos/packages/new_package that contains
# #np# a number of useful m4 files.  Any package that uses these m4 files
# #np# will have to add them to the repository.  The m4 files perform various
# #np# common checks and set common variable values to minimize redundant work.
#
AC_CONFIG_AUX_DIR(config)
#  #np# - Change file names in next line 
# Configure should create src/Teuchos_config.h from src/Teuchos_config.h.in

AM_CONFIG_HEADER(src/Teuchos_config.h:src/Teuchos_config.h.in)

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

# Specify required version of automake
AM_INIT_AUTOMAKE(1.10 no-define tar-ustar)

# Specify required version of autoconf.

AC_PREREQ(2.61)

# ------------------------------------------------------------------------
# Check if Fortran is enabled or not
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_OPTION(fortran, [enable Fortran support], FORTRAN_SUPPORT, yes)

# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI

#  #np# - can eliminate compiler checks below if your package does not use the
#         language corresponding to the check

# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ c++ cxx)
if test "X$ac_cv_use_fortran" != "Xno"; then
AC_PROG_F77(f77 g77 gfortran f90 xlf90 f95)
fi
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Check for perl
# ------------------------------------------------------------------------

# Check for perl to run scripts
AC_CHECK_PROG(HAVE_PERL, perl, yes, no)
AM_CONDITIONAL(USING_PERL, test X${HAVE_PERL} = Xyes)

#-------------------------------------------------------------------------
# See if the user is allowing us to assume that gnu-make will be used
#-------------------------------------------------------------------------

#TAC_ARG_ENABLE_FEATURE_SUB(
#  with, gnumake
#  ,[Specify that the the program GNUMake will be used to make all code.  If specified, perl must also be in your path!]
#  ,WITH_GNU_MAKE
#  ,no
#  )
#AM_CONDITIONAL(USING_WITH_GNU_MAKE, [test "X$ac_cv_use_with_gnumake" != "Xno"])

#AC_ARG_WITH(gnu-make,
# [AC_HELP_STRING([--with-gnu-make],
# [The program GNUMake will be used to make all code.  If specified, perl must also be in your path!])],
# [
#   USING_GNU_MAKE=1
# ],)

# ------------------------------------------------------------------------
# Checks for extra flags
# #np# - You can pull flag checks for languages not used in package 
# ------------------------------------------------------------------------
# Check if --with-flags present, prepend any specs to FLAGS

TAC_ARG_WITH_FLAGS(ccflags, CCFLAGS)
TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
TAC_ARG_WITH_FLAGS(cflags, CFLAGS)
if test "X$ac_cv_use_fortran" != "Xno"; then
TAC_ARG_WITH_FLAGS(fflags, FFLAGS)
fi
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

#  #np# - Customize the checks for headers needed in package being added.
# ------------------------------------------------------------------------
# Checks for header files - used in Teuchos_ConfigDefs.hpp
# ------------------------------------------------------------------------

AC_CHECK_HEADERS(
 [cerrno climits cstdlib cassert cmath iomanip iostream string  cstring complex algorithm numeric map list set vector sstream stdexcept cstdarg typeinfo memory utility],,
  AC_MSG_ERROR([[Cannot find standard header (see above).]])
)
# 2007/06/28: rabartl: Note that above you can not have a newline in the [...] 
# that gives the names of the header files!

AC_CHECK_HEADERS([fpu_control.h],,)

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE

# ------------------------------------------------------------------------
# Checks for compiler characteristics - Teuchos_ConfigDefs.h uses these, 
# so it makes sense to leave them in - see PORTABILITY in this directory
# ------------------------------------------------------------------------
AC_CXX_NAMESPACES
AC_CXX_HAVE_STL
AC_CXX_BOOL
AC_CXX_MUTABLE
AC_CXX_NUMERIC_LIMITS
AC_CXX_ABI_DEMANGLE
# ---------------------------------------------------------------------
# These are pretty obscure
# ---------------------------------------------------------------------
AC_CXX_NEW_FOR_SCOPING
AC_CXX_INVALID_TEMPLATE_QUALIFIER
AC_CXX_STD_IOS_BASE_FMTFLAGS
AC_CXX_STD_SPRINTF
AC_CXX_PROTECTED_NESTED_TEMPLATE_CLASS_ACCESS

# ------------------------------------------------------------------------
# Checks if tests should be built
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE(tests, [Make tests for all Trilinos packages buildable with 'make tests'], TESTS, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( teuchos, tests, [Make Teuchos tests buildable with 'make tests'], TEUCHOS_TESTS)
AM_CONDITIONAL(BUILD_TESTS, [test "X$ac_cv_use_teuchos_tests" != "Xno" ])

# ------------------------------------------------------------------------
# Checks if examples should be built
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE(examples, [Make examples for all Trilinos packages buildable with 'make examples'], EXAMPLES, yes)

TAC_ARG_ENABLE_FEATURE_SUB_CHECK( teuchos, examples, [Make Teuchos examples buildable with 'make examples'], TEUCHOS_EXAMPLES)
AM_CONDITIONAL(BUILD_EXAMPLES, [test "X$ac_cv_use_teuchos_examples" != "Xno"])

#We now build tests and examples through separate make targets, rather than
#during "make".  We still need to conditionally include the test and example
#in SUBDIRS, even though SUB_TEST and SUB_EXAMPLE will never be
#defined, so that the tests and examples are included in the distribution
#tarball.
AM_CONDITIONAL(SUB_TEST, test "X$ac_cv_use_sub_test" = "Xyes")
AM_CONDITIONAL(SUB_EXAMPLE, test "X$ac_cv_use_sub_example" = "Xyes")

TAC_ARG_ENABLE_FEATURE(libcheck, [Check for some third-party libraries including BLAS and LAPACK.  Cannot be disabled unless tests and examples are also disabled.  Cannot be disabled when fortran is disabled.], LIBCHECK, yes)

# ------------------------------------------------------------------------
# Specify other directories 
# ------------------------------------------------------------------------

# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS

# ------------------------------------------------------------------------
# Check for library specification errors
# -----------------------------------------------------------------------
AC_MSG_CHECKING(for library specification errors )
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_LINK_IFELSE([ AC_LANG_PROGRAM([[ ]], [[ ]]) ],
   [AC_MSG_RESULT(no errors detected)],
   [AC_MSG_RESULT(error found)
            echo "-----"
            echo "Cannot link with the current set of specified"
            echo "libraries and paths.  Look for typos and other"
            echo "possible errors in the values of --with-libs,"
            echo "--with-ldflags and other places where libraries"
            echo "are listed.  Grep for -"
            echo "checking for library specification errors"
            echo "-in packages/teuchos/config.log and look at"
            echo "the next several lines for a hint as to"
            echo "what the issue might be."
            echo "-----"
            AC_MSG_ERROR(cannot link)])
AC_LANG_RESTORE

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

# ------------------------------------------------------------------------
# Check to see if array boundschecking has been enabled
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE_SUB( teuchos, abc, [Enable array bounds checking], 
		TEUCHOS_ARRAY_BOUNDSCHECK, no)

#
#  --enable-gcc-checkedstl
#

TAC_ARG_ENABLE_FEATURE_SUB( gcc, checkedstl, [Enable the checked STL for compatible gcc compilers], GCC_CHECKEDSTL, no)
AM_CONDITIONAL(HAVE_GCC_CHECKEDSTL, [test "X$ac_cv_use_gcc_checkedstl" != "Xno"])

# ------------------------------------------------------------------------
# Checks for addons 
# ------------------------------------------------------------------------

#
#  --enable-teuchos-explicit_instantiation
#

TAC_ARG_ENABLE_FEATURE_SUB( teuchos, explicit_instantiation,
  [Enable explicit template instantiation for selected types],
   TEUCHOS_EXPLICIT_INSTANTIATION, no)
AM_CONDITIONAL(HAVE_TEUCHOS_EXPLICIT_INSTANTIATION,
  [test "X$ac_cv_use_teuchos_explicit_instantiation" != "Xno"])

#
#  --enable-teuchos-float
#

TAC_ARG_ENABLE_FEATURE_SUB( teuchos, float, [Enable Teuchos float code], TEUCHOS_FLOAT, yes)
AM_CONDITIONAL(HAVE_TEUCHOS_FLOAT, [test "X$ac_cv_use_teuchos_float" != "Xno"])

#
#  --enable-teuchos-complex
#

TAC_ARG_ENABLE_FEATURE_SUB( teuchos, complex, [Enable Teuchos complex code], TEUCHOS_COMPLEX, yes)
AM_CONDITIONAL(HAVE_TEUCHOS_COMPLEX, [test "X$ac_cv_use_teuchos_complex" != "Xno"])

#
#  --enable-teuchos-long-long
#

TAC_ARG_ENABLE_FEATURE_SUB( teuchos, long_long, [Enable long long int support in Teuchos], TEUCHOS_LONG_LONG_INT, no)
AM_CONDITIONAL(HAVE_TEUCHOS_LONG_LONG_INT, [test "X$ac_cv_use_teuchos_long_long_int" != "Xno"])

#
#  --enable-teuchos-extended
#

TAC_ARG_ENABLE_FEATURE_SUB( teuchos, extended, [Enable Teuchos extended code], TEUCHOS_EXTENDED, yes)
AM_CONDITIONAL(HAVE_TEUCHOS_EXTENDED, [test "X$ac_cv_use_teuchos_extended" != "Xno"])

#
#  --enable-teuchos-arprec 
# 
 
TAC_ARG_ENABLE_FEATURE_SUB( teuchos, arprec, [Enable Teuchos interface support for the ARbitrary PRECision Library], 
		TEUCHOS_ARPREC, no)  

#
#  --enable-teuchos-blasfloat
#

TAC_ARG_ENABLE_FEATURE_SUB( teuchos, blasfloat, [Enable support for float BLAS wrappers code], TEUCHOS_BLASFLOAT, yes)
AM_CONDITIONAL(HAVE_TEUCHOS_BLASFLOAT, [test "X$ac_cv_use_teuchos_blasfloat" != "Xno"])

if test "X$ac_cv_use_teuchos_arprec" != "Xno"; then

#
#  Check for ARPREC headers.
#
     AC_MSG_CHECKING(for ARPREC headers (mp/mpreal.h) ) 
     AC_TRY_CPP([#include "mp/mpreal.h"],
	[AC_MSG_RESULT(yes)],  
	[
	 AC_MSG_RESULT(no)
	 echo "-----"
     	 echo "Cannot find ARPREC headers."
     	 echo "Use --with-incdirs to specify ARPREC include directory."
     	 echo "-----"
     	 AC_MSG_ERROR(ARPREC headers cannot be found)
     	])
#
#  Check for ARPREC library.
#
     AC_MSG_CHECKING(for ARPREC library (libmp.a) )
     AC_LANG_SAVE
     AC_LANG_CPLUSPLUS
     AC_LINK_IFELSE([ AC_LANG_PROGRAM([[#include "mp/mpreal.h"]], 
	[[mp::mp_init(200);]]) ],
	[AC_MSG_RESULT(yes)],
	[AC_MSG_RESULT(no)
		 echo "-----"
	     	 echo "Cannot find ARPREC library."
     		 echo "Use --with-ldflags to specify ARPREC library directory."
	     	 echo "Use --with-libs to specify ARPREC library (-lmp)."
	     	 echo "-----"
	     	 AC_MSG_ERROR(ARPREC test cannot link)])
     AC_LANG_RESTORE

fi

#
#  --enable-teuchos-gmp
# 
 
TAC_ARG_ENABLE_FEATURE_SUB( teuchos, gmp, [Enable Teuchos interface support for the GNU MP Library], 
		TEUCHOS_GNU_MP, no)  

if test "X$ac_cv_use_teuchos_gmp" != "Xno"; then

#  TAC_ARG_WITH( trilinos3pldir, [Set Trilinos3PL directory], trilinos3pldir_is_not_specified )
#  if test "X${tac_with_teuchos_trilinos3pldir}" != "Xtrilinos3pldir_is_not_specified"; then
#
#  Check for GNU MP headers.
#
     AC_MSG_CHECKING(for GNU MP headers (gmpxx.h) ) 
     AC_TRY_CPP([#include "gmpxx.h"],
	[AC_MSG_RESULT(yes)],  
	[
	 AC_MSG_RESULT(no)
	 echo "-----"
     	 echo "Cannot find GNU MP headers."
     	 echo "Use --with-incdirs to specify GNU MP include directory."
     	 echo "-----"
     	 AC_MSG_ERROR(GNU MP headers cannot be found)
     	])
#
#  Check for GNU MP library.
#
     AC_MSG_CHECKING(for GNU MP library (libgmp.a / libgmpxx.a) )
     AC_LANG_SAVE
     AC_LANG_CPLUSPLUS
     AC_LINK_IFELSE([ AC_LANG_PROGRAM([[#include "gmpxx.h" 
					#include "gmp.h"]], 
					[[mpf_class one = 1.0; sqrt(one);]]) ],
	[AC_MSG_RESULT(yes)],
	[AC_MSG_RESULT(no)
		 echo "-----"
	     	 echo "Cannot find GNU MP library."
     		 echo "Use --with-ldflags to specify GNU MP library directory."
	     	 echo "Use --with-libs to specify GNU MP library (-lgmp -lgmpxx)."
	     	 echo "-----"
	     	 AC_MSG_ERROR(GNU MP test cannot link)])
     AC_LANG_RESTORE

fi

#
#  --enable-teuchos-qd
# 
 
TAC_ARG_ENABLE_FEATURE_SUB( teuchos, qd, [Enable Teuchos interface support for the QD library], 
		TEUCHOS_QD, no)  

if test "X$ac_cv_use_teuchos_qd" != "Xno"; then

  TAC_ARG_ENABLE_FEATURE_SUB( teuchos, doubletodd, [Instruct Teuchos ScalarTraits to promote double to dd_real], TEUCHOS_DOUBLE_TO_DD, no)

#  TAC_ARG_WITH( trilinos3pldir, [Set Trilinos3PL directory], trilinos3pldir_is_not_specified )
#  if test "X${tac_with_teuchos_trilinos3pldir}" != "Xtrilinos3pldir_is_not_specified"; then
#
#  Check for QD headers.
#
     AC_MSG_CHECKING(for QD headers (qd/qd_real.h) ) 
     AC_TRY_CPP([#include "qd/qd_real.h"],
	[AC_MSG_RESULT(yes)],  
	[
	 AC_MSG_RESULT(no)
	 echo "-----"
     	 echo "Cannot find QD headers."
     	 echo "Use --with-incdirs to specify QD include directory."
     	 echo "-----"
     	 AC_MSG_ERROR(QD headers cannot be found)
     	])
#
#  Check for QD library.
#
     AC_MSG_CHECKING(for QD library (libqd.a) )
     AC_LANG_SAVE
     AC_LANG_CPLUSPLUS
     AC_LINK_IFELSE([ AC_LANG_PROGRAM([[#include "qd/qd_real.h"]], 
					[[qd_real one(1.0); sqrt(one);]]) ],
	[AC_MSG_RESULT(yes)],
	[AC_MSG_RESULT(no)
		 echo "-----"
	     	 echo "Cannot find QD library."
     		 echo "Use --with-ldflags to specify QD library directory."
	     	 echo "Use --with-libs to specify QD library (-lqd)."
	     	 echo "-----"
	     	 AC_MSG_ERROR(QD test cannot link)])
     AC_LANG_RESTORE

fi

#
#  --enable-teuchos-expat 
# 
 
TAC_ARG_ENABLE_FEATURE_SUB( teuchos, expat, [Enable Teuchos interface support for the Expat XML Library], 
		TEUCHOS_EXPAT, no)  

if test "X$ac_cv_use_teuchos_expat" != "Xno"; then

#  TAC_ARG_WITH( trilinos3pldir, [Set Trilinos3PL directory], trilinos3pldir_is_not_specified )
#  if test "X${tac_with_teuchos_trilinos3pldir}" != "Xtrilinos3pldir_is_not_specified"; then
#
#  Check for Expat headers.
#
     AC_MSG_CHECKING(for Expat headers (expat.h) ) 
     AC_TRY_CPP([#include "expat.h"],
	[AC_MSG_RESULT(yes)],  
	[
	 AC_MSG_RESULT(no)
	 echo "-----"
     	 echo "Cannot find Expat headers."
     	 echo "Use --with-incdirs to specify Expat include directory."
     	 echo "-----"
     	 AC_MSG_ERROR(Expat headers cannot be found)
     	])
#
#  Check for Expat library.
#
     AC_MSG_CHECKING(for Expat library (libexpat.a) )
     AC_LANG_SAVE
     AC_LANG_CPLUSPLUS
     AC_LINK_IFELSE([ AC_LANG_PROGRAM([[#include "expat.h"]], [[XML_Parser parser = XML_ParserCreate(NULL);]]) ],
	[AC_MSG_RESULT(yes)],
	[LIBS="-lexpat $LIBS"
	 AC_LINK_IFELSE([ AC_LANG_PROGRAM([[#include "expat.h"]], [[XML_Parser parser = XML_ParserCreate(NULL);]]) ],
        [AC_MSG_RESULT(yes)],
        [AC_MSG_RESULT(no)
                 echo "-----"
                 echo "Cannot find Expat library."
                 echo "Use --with-ldflags to specify Expat library directory."
                 echo "Use --with-libs to specify Expat library (-lexpat)."
                 echo "-----"
                 AC_MSG_ERROR(Expat test cannot link)])
	])
     AC_LANG_RESTORE

fi

#
#  --enable-teuchos-libxml2 
# 
 
TAC_ARG_ENABLE_FEATURE_SUB( teuchos, libxml2, [Enable Teuchos interface support for the LibXML2 Library], 
		TEUCHOS_LIBXML2, no)  

if test "X$ac_cv_use_teuchos_libxml2" != "Xno"; then

#  TAC_ARG_WITH( trilinos3pldir, [Set Trilinos3PL directory], trilinos3pldir_is_not_specified )
#  if test "X${tac_with_teuchos_trilinos3pldir}" != "Xtrilinos3pldir_is_not_specified"; then
#
#  Check for LibXML2 headers.
#
     AC_MSG_CHECKING(for LibXML2 headers (libxml/parser.h) ) 
     CPPFLAGS="-I/usr/include/libxml2 $CPPFLAGS";
     AC_TRY_CPP([#include "libxml/parser.h"],
	[AC_MSG_RESULT(yes)],  
	[
	 AC_MSG_RESULT(no)
	 echo "-----"
     	 echo "Cannot find LibXML2 headers."
     	 echo "Use --with-incdirs to specify LibXML2 include directory."
     	 echo "-----"
     	 AC_MSG_ERROR(LibXML2 headers cannot be found)
     	])
#
#  Check for LibXML2 library.
#
     AC_MSG_CHECKING(for LibXML2 library (libxml2.a) )
     AC_LANG_SAVE
     AC_LANG_CPLUSPLUS
     AC_LINK_IFELSE([ AC_LANG_PROGRAM([[#include "libxml/tree.h"]], 
        [[xmlDocPtr doc = NULL; doc = xmlNewDoc(BAD_CAST "1.0"); xmlFreeDoc(doc);]]) ],
	[AC_MSG_RESULT(yes)],
	[LIBS="-lxml2 $LIBS"
	 AC_LINK_IFELSE([ AC_LANG_PROGRAM([[#include "libxml/tree.h"]], 
         [[xmlDocPtr doc = NULL; doc = xmlNewDoc(BAD_CAST "1.0"); xmlFreeDoc(doc);]]) ],
	[AC_MSG_RESULT(yes)],
	[AC_MSG_RESULT(no)
		 echo "-----"
	     	 echo "Cannot find LibXML2 library."
     		 echo "Use --with-ldflags to specify LibXML2 library directory."
	     	 echo "Use --with-libs to specify LibXML2 library (-lxml2)."
	     	 echo "-----"
	     	 AC_MSG_ERROR(LibXML2 test cannot link)])
	])
     AC_LANG_RESTORE

fi

if test "X${HAVE_TEUCHOS_LIBXML2}" = "Xyes"; then
  if test "X${HAVE_TEUCHOS_EXPAT}" = "Xyes"; then
    AC_MSG_ERROR([[You have specified both Expat and LibXML2. Please enable only one XML parser!]])
  fi  
fi


#
#  --enable-teuchos-boost
#

TAC_ARG_ENABLE_FEATURE_SUB( teuchos, boost, [Enable support for some of boost], TEUCHOS_BOOST, no)
AM_CONDITIONAL(HAVE_TEUCHOS_BOOST, [test "X$ac_cv_use_teuchos_boost" != "Xno"])

#
# --enable-teuchos-comm_timers
#

TAC_ARG_ENABLE_FEATURE_SUB( teuchos, comm_timers, [Enable timers for Teuchos::Comm related software], 
		TEUCHOS_COMM_TIMERS, yes )
#AM_CONDITIONAL(HAVE_TEUCHOS_COMM_TIMERS, [test "X$ac_cv_use_teuchos_comm_timers" != "Xno"])

#
# --enable-teuchos-debug
#

TAC_ARG_ENABLE_FEATURE_SUB( teuchos, debug, [Enable optional extra debugging checks], 
		TEUCHOS_DEBUG, no )

#
# --enable-teuchos-demangle
#

TAC_ARG_ENABLE_FEATURE_SUB( teuchos, demangle,
  [Name demangling for types (g++ only).  By default, configure will try to figure out if this is supported automatically but this does not work on every platoform.  Use --disable-teuchos-demangle to turn this off and avoid problems],
  TEUCHOS_DEMANGLE, yes )

# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

# Define F77_FUNC that will be used to link with Fortran subroutines.
if test "X$ac_cv_use_fortran" != "Xno"; then
AC_F77_WRAPPERS
fi

# If tests, examples and libcheck are disabled, we don't have to check
# for these libraries.
if test "X$ac_cv_use_teuchos_tests" != "Xno" || test "X$ac_cv_use_teuchos_examples" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then

# #np# - 
# To add any libraries that you want to check for, replace pow below 
# with the name of an entry point in that library.
# checks for libraries now
dnl Replace `main' with a function in -lm:
#  #np# These lines can be removed if the math lib is not needed.
AC_SEARCH_LIBS(pow,[m],,[AC_MSG_ERROR(Error encountered when checking for math library.  ***Almost 
always this error is caused by a typo in the list of configure arguments 
that is not related to the math library.  Look especially for missing 
quotes or hyphens in any value related to compiler flags or libraries.
For an error message that may assist in finding the problem grep for 
"checking for library containing pow" in the file 
packages/teuchos/config.log and read the next several lines.)])
AC_CHECK_FUNCS([pow sqrt])

### Include Blas/Lapack example here. ###  For now, developers should see
### The configure.ac file for epetra (Trilinos/packages/epetra) if blas and
### or Lapack are needed.

TAC_BLAS([], [AC_MSG_ERROR([[Cannot find blas library, specify a path using --with-blas=DIR/LIB (for example --with-blas=/usr/path/lib/libcxml.a)]])])
TAC_LAPACK([], [AC_MSG_ERROR([[Cannot find lapack library, specify a path using --with-lapack=DIR/LIB (ex. --with-lapack=/usr/path/lib/libcxml.a)]])])

# Find valid LAPACK _GGEV function [LAPACK library may only offer deprecated _GEGV function].
ACX_GGEV

else

  if test "X$ac_cv_use_fortran" != "Xno"; then
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

fi

# check whether C/C++ has problems calling cdotc and zdotc and whether we know how to fix it
TAC_ARG_ENABLE_FEATURE_SUB( teuchos, dotcheck, [Enable checking on complex BLAS dot routines], 
		TEUCHOS_DOTCHECK, no )
if test "X$ac_cv_use_teuchos_dotcheck" != "Xno"; then
        AC_CXX_COMPLEX_BLAS_PROBLEM
        if test "X$ac_cv_cxx_complex_blas_problem" == "Xyes"; then
                AC_CXX_COMPLEX_BLAS_PROBLEM_CAN_BE_FIXED
        fi
fi

# end of the list of libraries that don't need to be checked for if
# tests and examples are disabled.

# ------------------------------------------------------------------------
# Checks for Makefile.export related systems
# ------------------------------------------------------------------------
TAC_ARG_ENABLE_EXPORT_MAKEFILES(yes)

# ------------------------------------------------------------------------
# Checks for linker characteristics
# ------------------------------------------------------------------------

# Determine libraries needed for linking with Fortran
#if test "X$ac_cv_use_fortran" != "Xno"; then
AC_F77_LIBRARY_LDFLAGS
#fi

#  ----------------------------------------------------
#  Needs this because Dec handles templates differently
#  ----------------------------------------------------
#  #np# - Only uncomment the following lines if the package
#         creates a subdirectory called cxx_repository for
#         template object files on DEC.  (HINT: try to build on 
#         DEC/CPLANT and see if template-related errors occur)

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
# Get the absolute path to the build directory.  
# 
# We need this as a reference point from which 
# to figure out paths to the libraries for the other Trilinos packages. 
# 
# ------------------------------------------------------------------------

AC_MSG_CHECKING([finding build directory])
if test "X" = "X"; then
  TEUCHOS_BUILD_DIRECTORY=`pwd`
  AC_MSG_RESULT($TEUCHOS_BUILD_DIRECTORY)
fi


# ------------------------------------------------------------------------
# Perform substitutions in output files
# ------------------------------------------------------------------------

AC_SUBST(ac_aux_dir)
AC_SUBST(TEUCHOS_BUILD_DIRECTORY)
AC_SUBST(MPI_LIBS)
AC_SUBST(MPI_LIBDIR)
AC_SUBST(MPI_INC)
AC_SUBST(MPI_DIR)
AC_SUBST(AMDEP_TRUE)
AC_SUBST(AMDEP_FALSE)


# ------------------------------------------------------------------------
# Output files
# ------------------------------------------------------------------------
##
#  #np# - Most packages should have an example and a full suite 
#  of tests.  New package comes only with src and example.
#  You will need to change AC_CONFIG_FILES below and Makefile.am
#  to add a new directory.
AC_CONFIG_FILES([
		Makefile 
		Makefile.export 
		Makefile.export.teuchos 
		src/Makefile
		test/BLAS/Makefile
		test/Comm/Makefile
		test/Containers/Makefile
		test/DenseMatrix/Makefile
		test/FancyOutputting/Makefile
		test/Handle/Makefile
		test/LAPACK/Makefile
		test/Makefile
		test/MemoryManagement/Makefile
		test/ObjectBuilder/Makefile
		test/ParameterList/Makefile
		test/Polling/Makefile
		test/Polynomial/Makefile
		test/ScalarTraits/Makefile
		test/Time/Makefile
		test/UnitTest/Makefile
		test/Workspace/Makefile
		test/XML/Makefile
		test/XMLParser/Makefile
		test/dyn_cast/Makefile
		test/exit/Makefile
		test/scripts/Makefile
		test/scripts/daily/Makefile
		test/scripts/daily/mpi/Makefile
		test/scripts/daily/serial/run-teuchos-tests.pl
		example/Makefile
		example/BLAS/Makefile
		example/CommandLineProcessor/Makefile
		example/DenseMatrix/Makefile
		example/hilbert/Makefile
		example/LAPACK/Makefile
		example/ParameterList/Makefile
		example/RefCountPtr/Makefile
		example/TimeMonitor/Makefile
		example/XML/Makefile
		])

AC_OUTPUT()

# Bye World!
#  #np# change package name
echo "-----------------------------------------"
echo "Finished Running Teuchos Configure Script"
echo "-----------------------------------------"
