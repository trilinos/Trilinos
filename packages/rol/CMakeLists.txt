INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)


#
# A) Define the package
#

TRIBITS_PACKAGE(ROL)

#
# B) Set up package-specific options
#

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG
  HAVE_ROL_DEBUG
  "Enable a host of runtime debug checking."
  ${${PROJECT_NAME}_ENABLE_DEBUG}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_TIMERS
  ROL_TIMERS
  "Build ROL with Teuchos TimeMonitors enabled."
  OFF
  )

# Build Options
SET( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
INCLUDE(BuildOptions)

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

IF( ROL_ENABLE_Sacado ) 
ADD_SUBDIRECTORY(adapters/sacado)
ENDIF()

IF (ROL_ENABLE_ArrayFireCPU)
  ADD_SUBDIRECTORY(adapters/arrayfire)
ENDIF()

# Adapters which require use of Teuchos

GET_PROPERTY( PTR_STRING GLOBAL PROPERTY PTR_IMPL )
GET_PROPERTY( PARAMETERLIST_STRING GLOBAL PROPERTY PARAMETERLIST_IMPL )
 
IF( PTR_STRING STREQUAL "Teuchos::RCP" )
IF( PARAMETERLIST_STRING STREQUAL "Teuchos::ParameterList" )

ADD_SUBDIRECTORY(adapters/teuchos)

IF( ROL_ENABLE_Epetra )
ADD_SUBDIRECTORY(adapters/epetra)
ENDIF()

IF( ROL_ENABLE_Tpetra )
ADD_SUBDIRECTORY(adapters/tpetra)
ENDIF()

IF( ROL_ENABLE_Belos )
ADD_SUBDIRECTORY(adapters/belos)
ENDIF()

IF( ROL_ENABLE_Thyra )
ADD_SUBDIRECTORY(adapters/thyra)
ENDIF()

IF( ROL_ENABLE_TriKota )
ADD_SUBDIRECTORY(adapters/trikota)
ENDIF()

IF( ROL_ENABLE_PEBBL )
ADD_SUBDIRECTORY(adapters/pebbl)
ENDIF()

IF (ROL_ENABLE_Boost AND ROL_ENABLE_MiniTensor)
  ADD_SUBDIRECTORY(adapters/minitensor)
ENDIF()


ENDIF( PARAMETERLIST_STRING STREQUAL "Teuchos::ParameterList" )
ENDIF( PTR_STRING           STREQUAL "Teuchos::RCP" )

IF (ROL_ENABLE_Eigen)
  ADD_SUBDIRECTORY(adapters/eigen)
ENDIF()


TRIBITS_ADD_TEST_DIRECTORIES(test)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(tutorial)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
