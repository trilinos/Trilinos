INCLUDE(PackageMacros)

#
# A) Define the package
#

PACKAGE(Amesos2)

SET(${PACKAGE_NAME_UC}_VERSION "0.3d")
SET(${PACKAGE_NAME_UC}_RELEASE_DATE "07/28/2011")

#
# B) Set up package-specific options
#

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_BUILD_EXPERIMENTAL
  HAVE_EXPERIMENTAL
  "Enable experimental features"
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_VERBOSE_DEBUG
  HAVE_${PACKAGE_NAME_UC}_VERBOSE_DEBUG
  "Enable verbose output during debug-mode"
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_TIMERS
  HAVE_${PACKAGE_NAME_UC}_TIMERS
  "Enable performance timers"
  OFF )

#
# TPL Solver interfaces are enabled automatically based on Trilinos
# TPL configuration.
#
# Config options such as Amesos2_ENABLE_XXXXXX are automatically
# generated in case a user would wish to explicitely disable Amesos2
# support for an otherwise enabled third-part solver library
#

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_LAPACK
  HAVE_${PACKAGE_NAME_UC}_LAPACK
  "Enable the dense serial solver interface"
  OFF )

PACKAGE_ADD_DEBUG_OPTION()

PACKAGE_ADD_EXPLICIT_INSTANTIATION_OPTION()

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
