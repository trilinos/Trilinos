# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# @HEADER

# ------------------------------------------------------------------------
# Initialization 
# ------------------------------------------------------------------------



# This must be the first line in configure.ac.
# Optional 3rd argument is email address for bugs.

#  #np# - package name, version number, and e-mail address below
AC_INIT(sample, 1.0d, jmwille@sandia.gov)

echo "----------------------------------------"
echo "Running sample Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.  Any file in that directory will do,
# preferably one that is unlikely to be removed or renamed.

# #auto np# - file name below
#AC_CONFIG_SRCDIR([src/Newp_Hello.cpp])

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

#
# #np# This a a subdirectory of Trilinos/packages/new_package that contains
# #np# a number of useful m4 files.  Any package that uses these m4 files
# #np# will have to add them to the repository.  The m4 files perform various
# #np# common checks and set common variable values to minimize redundant work.
#
AC_CONFIG_AUX_DIR(config)
#  #auto np# - Change file names in next line 
# Configure should create src/New_Package_config.h from src/New_Package_config.h.in

#AM_CONFIG_HEADER(src/New_Package_config.h:src/New_Package_config.h.in)

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

#  - Required version of automake.
AM_INIT_AUTOMAKE(1.10 no-define tar-ustar)

# Specify required version of autoconf.

AC_PREREQ(2.61)

# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

#AC_PROG_CC(cc gcc)
#AC_PROG_CXX(CC g++ gcc c++ cxx)
#AC_PROG_F77(f77 g77 f90 xlf90)
#AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Checks for extra flags
# #np# - You can pull flag checks for languages not used in package.  Please 
#	  note that if you use F77_FUNC to determine Fortran name mangling, 
#	  you should not remove the Fortran compiler check or the check for 
#	  Fortran flags.  Doing so will prevent the detection of the proper 
#	  name mangling in some cases.
# ------------------------------------------------------------------------
# Check if --with-flags present, prepend any specs to FLAGS

#TAC_ARG_WITH_FLAGS(ccflags, CCFLAGS)
#TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
#TAC_ARG_WITH_FLAGS(cflags, CFLAGS)
#TAC_ARG_WITH_FLAGS(fflags, FFLAGS)
#TAC_ARG_WITH_LIBS
#TAC_ARG_WITH_FLAGS(ldflags, LDFLAGS)

# ------------------------------------------------------------------------
# Alternate archiver
# ------------------------------------------------------------------------

#TAC_ARG_WITH_AR

# ------------------------------------------------------------------------
# Checks for header files - used in Epetra_ConfigDefs.h 
# ------------------------------------------------------------------------

#AC_CHECK_HEADERS([cstdlib stdlib.h], [break],)

#if test "X$ac_cv_header_cstdlib" != "Xyes" && test "X$ac_cv_header_stdlib_h" != "Xyes"; then
#  AC_MSG_ERROR([[Cannot find cstdlib or stdlib.h.  New_Package requires at least one or the other.]])
#fi

# Checks for typedefs, structures, and compiler characteristics.
#AC_C_INLINE

# ------------------------------------------------------------------------
# Checks for compiler characteristics - Epetra_ConfigDefs.h uses these, 
# so it makes sense to leave them in - see PORTABILITY in this directory
# ------------------------------------------------------------------------
#AC_CXX_NAMESPACES
#AC_CXX_HAVE_STL
#AC_CXX_BOOL
#AC_CXX_MUTABLE
# ---------------------------------------------------------------------
# These are pretty obscure
# ---------------------------------------------------------------------
#AC_CXX_NEW_FOR_SCOPING
#AC_CXX_STD_IOS_BASE_FMTFLAGS
#AC_CXX_STD_SPRINTF

# ------------------------------------------------------------------------
# Checks for Makefile.export related systems
# ------------------------------------------------------------------------
#TAC_ARG_ENABLE_EXPORT_MAKEFILES(yes)

# ------------------------------------------------------------------------
# Checks if tests and examples should be built
# ------------------------------------------------------------------------

#  #np# - These options can disable the tests and examples of a package.
#  #np# - Packages that do not have tests or examples should #-out the 
#  #np# - option(s) that does (do) not apply.

#TAC_ARG_ENABLE_FEATURE(tests, [Make tests for all Trilinos packages buildable with 'make tests'], TESTS, yes)
#TAC_ARG_ENABLE_FEATURE_SUB_CHECK( new_package, tests, [Make New_Package tests buildable with 'make tests'], NEW_PACKAGE_TESTS)
#AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_new_package_tests" != "Xno")

#TAC_ARG_ENABLE_FEATURE(examples, [Make examples for all Trilinos packages buildable with 'make examples'], EXAMPLES, yes)
#TAC_ARG_ENABLE_FEATURE_SUB_CHECK( new_package, examples, [Make New_Package examples buildable with 'make examples'], NEW_PACKAGE_EXAMPLES)

#We now build tests and examples through separate make targets, rather than
#during "make".  We still need to conditionally include the test and example
#in SUBDIRS, even though SUB_TEST and SUB_EXAMPLE will never be
#defined, so that the tests and examples are included in the distribution
#tarball.
#AM_CONDITIONAL(SUB_TEST, test "X$ac_cv_use_sub_test" = "Xyes")
#AM_CONDITIONAL(SUB_EXAMPLE, test "X$ac_cv_use_sub_example" = "Xyes")

# #np# - Yet another opportunity to remove code if you aren't
# using Fortran
# Define F77_FUNC that will be used to link with Fortran subroutines.
#AC_F77_WRAPPERS

# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

#  #np# These lines can be removed if the math lib is not needed.
#AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
#AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))

#fi
# end of the list of libraries that don't need to be checked for if
# tests and examples are disabled.

# ------------------------------------------------------------------------
# Checks for linker characteristics
# ------------------------------------------------------------------------

# Determine libraries needed for linking with Fortran
#AC_F77_LIBRARY_LDFLAGS

# ------------------------------------------------------------------------
# Perform substitutions in output files
# ------------------------------------------------------------------------

AC_SUBST(ac_aux_dir)

# ------------------------------------------------------------------------
# Output files
# ------------------------------------------------------------------------

AC_CONFIG_FILES([Makefile])

AC_OUTPUT()

# Bye World!
echo "----------------------------------------"
echo "Finished Running sample Configure Script"
echo "----------------------------------------"
