TRIBITS_SUBPACKAGE(DiscFE)

# Optional CAMAL TPL dependency
###############################
SET(PANZER_HAVE_CAMAL ${${PARENT_PACKAGE_NAME}_ENABLE_CAMAL} )

#Optional PAPI TPL dependency
#############################
COMBINED_OPTION( ${PACKAGE_NAME}_BUILD_PAPI_SUPPORT
  DEP_OPTIONS_NAMES ${PACKAGE_NAME}_ENABLE_PAPI Teuchos_ENABLE_LONG_LONG_INT
  DOCSTR "Enables support for the PAPI performance counters."
    " Defaults to ON if both the PAPI TPL is enabled and long long int "
    " support is enabled in Teuchos, otherwise defaults to OFF."
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

TRIBITS_SUBPACKAGE_POSTPROCESS()
