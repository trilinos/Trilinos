TRIBITS_SUBPACKAGE(DiscFE)

# Optional CAMAL TPL dependency
###############################
SET(PANZER_HAVE_CAMAL ${${PARENT_PACKAGE_NAME}_ENABLE_CAMAL} )

#Optional PAPI TPL dependency
#############################
COMBINED_OPTION( ${PACKAGE_NAME}_BUILD_PAPI_SUPPORT
  DEP_OPTIONS_NAMES ${PACKAGE_NAME}_ENABLE_PAPI
  DOCSTR "Enables support for the PAPI performance counters."
    " Defaults to ON if the PAPI TPL is enabled,"
    " otherwise defaults to OFF."
  )

#Optional Hessian support
#############################

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PARENT_PACKAGE_NAME}_ENABLE_HESSIAN_SUPPORT
  ${PARENT_PACKAGE_NAME}_BUILD_HESSIAN_SUPPORT
  "Enable building of Hessian AD types"
  OFF
  )

IF(${PARENT_PACKAGE_NAME}_BUILD_HESSIAN_SUPPORT)
   MESSAGE("-- Hessian support On")
ELSE()
   MESSAGE("-- Hessian support Off")
ENDIF()

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

TRIBITS_SUBPACKAGE_POSTPROCESS()
