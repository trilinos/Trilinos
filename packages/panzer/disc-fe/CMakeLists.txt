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


ADD_SUBDIRECTORY(src)
TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_SUBPACKAGE_POSTPROCESS()
