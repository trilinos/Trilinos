INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(CombinedOption)

#
# A) Define the package
#

TRIBITS_PACKAGE(Panzer)

#
# B) Set up package-specific options
#

TRIBITS_ADD_ENABLE_TEUCHOS_TIME_MONITOR_OPTION()

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG
  PANZER_DEBUG
  "Build with debug code enabled."
  OFF )

# optional package specific configurations
############################################################

IF(${PACKAGE_NAME}_ENABLE_SEACASIoss)
   GLOBAL_SET(HAVE_IOSS ON)
ELSE()
   GLOBAL_SET(HAVE_IOSS OFF)
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_Stokhos)
   GLOBAL_SET(HAVE_STOKHOS ON)
ELSE()
   GLOBAL_SET(HAVE_STOKHOS OFF)
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_Teko)
   GLOBAL_SET(HAVE_TEKO ON)
ELSE()
   GLOBAL_SET(HAVE_TEKO OFF)
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_MueLu)
   GLOBAL_SET(HAVE_MUELU ON)
ELSE()
   GLOBAL_SET(HAVE_MUELU OFF)
ENDIF()

############################################################

COMBINED_OPTION( ${PACKAGE_NAME}_BUILD_SOLVER_DRIVERS_EPETRA
  DEP_OPTIONS_NAMES ${PACKAGE_NAME}_ENABLE_Stratimikos ${PACKAGE_NAME}_ENABLE_Piro ${PACKAGE_NAME}_ENABLE_NOX ${PACKAGE_NAME}_ENABLE_Rythmos
  DOCSTR "Enables support for the solver drivers for the tests/examples."
    " Defaults to ON if Stratimikos, Piro, NOX, and Rythmos are"
    " enabled, otherwise defaults to OFF."
  )

COMBINED_OPTION( ${PACKAGE_NAME}_BUILD_STK_AND_SEACAS_SUPPORT
  DEP_OPTIONS_NAMES ${PACKAGE_NAME}_ENABLE_STK ${PACKAGE_NAME}_ENABLE_SEACASIoss ${PACKAGE_NAME}_ENABLE_SEACASExodus
  DOCSTR "Enables support for the STK mesh database and SEACAS I/O Layers."
    " Defaults to ON if STK, SEACASIoss, SEACASExodus are"
    " enabled, otherwise defaults to OFF."
  )

COMBINED_OPTION( ${PACKAGE_NAME}_BUILD_PAPI_SUPPORT
  DEP_OPTIONS_NAMES ${PACKAGE_NAME}_ENABLE_PAPI Teuchos_ENABLE_LONG_LONG_INT
  DOCSTR "Enables support for the PAPI performance counters."
    " Defaults to ON if both the PAPI TPL is enabled and long long int "
    " support is enabled in Teuchos, otherwise defaults to OFF."
  )

# Explicit Template Instantiation
TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()
SET(${PARENT_PACKAGE_NAME}_TEMPLATE_DECL_EXT "_decl.hpp")
SET(${PARENT_PACKAGE_NAME}_TEMPLATE_DEF_EXT "_impl.hpp")

## One day hope to add enabling evaluation types from configure flag
#ADD_OPTION(${PACKAGE_NAME}_EXPLICIT_INSTANTIATION_EVALUATION_TYPES
#  "Set the explicit instantiation evaluation types."
#  "panzer::Traits::Residual 2 panzer::Traits::Jacobian 1" )
#FILE(WRITE ${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE} ${DEPS_XML} )


#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(adapters)

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
