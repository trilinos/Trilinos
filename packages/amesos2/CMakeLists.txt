INCLUDE(PackageMacros)

#
# A) Define the package
#

PACKAGE(Amesos2)

SET(${PACKAGE_NAME_UC}_VERSION "0.1d")
SET(${PACKAGE_NAME_UC}_RELEASE_DATE "06/06/2011")

#
# B) Set up package-specific options
#

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_BUILD_EXPERIMENTAL
  HAVE_EXPERIMENTAL
  "Enable experimental features"
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_VERBOSE_DEBUG
  HAVE_${PACKAGE_NAME_UC}_VERBOSE_DEBUG
  "Enable output during debug-mode"
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Timers
  HAVE_${PACKAGE_NAME_UC}_TIMERS
  "Enable performance timers"
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_EpetraExt
  HAVE_${PACKAGE_NAME_UC}_EPETRAEXT
  "Enable EpetraExt functionality (for examples and testing)"
# ${Trilinos_ENABLE_EpetraExt} )  # How to set to ON if enabled in Trilinos?
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Epetra
  HAVE_${PACKAGE_NAME_UC}_EPETRA
  "Enable support for Epetra matrix and multivector objects"
#  ${Trilinos_ENABLE_Epetra} )
  OFF )


#
# Add options to enable TPL solver interfaces
#

# ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_KLU2
#   HAVE_AMESOS2_KLU2
#   "Enable KLU2 functionality."
#   ON )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_SuperLU
  HAVE_${PACKAGE_NAME_UC}_SUPERLU
  "Enable SuperLU functionality.  Requires the SuperLU TPL."
  ON )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_SuperLUMT
  HAVE_${PACKAGE_NAME_UC}_SUPERLUMT
  "Enable SuperLU_MT functionality.  Requires the SuperLUMT TPL"
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_SuperLUDist
  HAVE_${PACKAGE_NAME_UC}_SUPERLUDIST
  "Enable SuperLU_DIST functionality.  Requires the SuperLUDist TPL"
  OFF )


PACKAGE_ADD_DEBUG_OPTION()

# PACKAGE_ADD_EXPLICIT_INSTANTIATION_OPTION()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
