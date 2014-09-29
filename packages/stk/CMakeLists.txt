
#
# A) Define the package
#

TRIBITS_PACKAGE_DECL( STK  )

TRIBITS_ADD_DEBUG_OPTION()
TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

#
# B) Set up package-specific options
#

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_Trilinos
  HAVE_STK_Trilinos
  "Building as a Trilinos package"
  ON
  )

# set fortran mangling options

IF("${FC_FN_UNDERSCORE}" STREQUAL "NO_UNDER")
  SET(FORTRAN_NO_UNDERSCORE ON)
ELSEIF("${FC_FN_UNDERSCORE}" STREQUAL "UNDER")
  SET(FORTRAN_ONE_UNDERSCORE ON)
ELSEIF("${FC_FN_UNDERSCORE}" STREQUAL "SECOND_UNDER")
  SET(FORTRAN_TWO_UNDERSCORES ON)
ELSE()
  MESSAGE("Could not determine the Fortran mangling defaulting to no underscore.")
  SET(FORTRAN_NO_UNDERSCORE ON)
ENDIF()


#
# C) Add the libraries, tests, and examples
TRIBITS_PROCESS_SUBPACKAGES()

TRIBITS_PACKAGE_DEF()

#need to find out what to exclude
#TRIBITS_EXCLUDE_FILES(
#  )
  
#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()

