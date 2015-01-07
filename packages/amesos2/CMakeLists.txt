
#
# A) Define the package
#

TRIBITS_PACKAGE(Amesos2)

SET(${PACKAGE_NAME_UC}_VERSION "0.3d")
SET(${PACKAGE_NAME_UC}_RELEASE_DATE "07/28/2011")

#
# B) Set up package-specific options
#


TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_VERBOSE_DEBUG
  HAVE_${PACKAGE_NAME_UC}_VERBOSE_DEBUG
  "Enable verbose output during debug-mode"
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_TIMERS
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


IF(${PACKAGE_NAME}_ENABLE_Cholmod)
  SET(HAVE_AMESOS2_CHOLMOD ON )
ENDIF()


TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_LAPACK
  HAVE_${PACKAGE_NAME_UC}_LAPACK
  "Enable the dense serial solver interface"
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_KLU2
  HAVE_AMESOS2_KLU2
  "Enable KLU2 in Amesos2"
  OFF
)

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Basker
   HAVE_AMESOS2_BASKER
   "Enable Basker in Amesos2"
   ON
)

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_Experimental
  HAVE_AMESOS2_EXPERIMENTAL
  "Enable building and installation of experimental Amesos2 features."
  NO )

TRIBITS_ADD_DEBUG_OPTION()

TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

#leaving in place for future excludes
#TRIBITS_EXCLUDE_FILES(
#  )

TRIBITS_PACKAGE_POSTPROCESS()
