#
# A) Define the package
#

TRIBITS_PACKAGE(Stokhos)

#
# B) Set up package-specific options
#

TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG
  STOKHOS_DEBUG
  "Enable debug code in stokhos"
  ${${PROJECT_NAME}_ENABLE_DEBUG} )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_TEUCHOS_TIME_MONITOR
  STOKHOS_TEUCHOS_TIME_MONITOR
  "Build with Teuchos TimeMonitors enabled."
  "${${PROJECT_NAME}_ENABLE_TEUCHOS_TIME_MONITOR}" )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Experimental
  STOKHOS_EXPERIMENTAL
  "Enable experimental code in stokhos"
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Ensemble_Reduction
  HAVE_STOKHOS_ENSEMBLE_REDUCT
  "Enable reduction across across ensemble components."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Stokhos_ENABLE_CXX11
  HAVE_STOKHOS_CXX11
  "Enable C++11 support in Stokhos."
  "${${PROJECT_NAME}_ENABLE_CXX11}"
  )

SET(Stokhos_ENABLE_Ensemble_Scalar_Type_Default OFF)
IF(Stokhos_ENABLE_Sacado)
  SET(Stokhos_ENABLE_Ensemble_Scalar_Type_Default ON)
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(Stokhos_ENABLE_Ensemble_Scalar_Type
  HAVE_STOKHOS_ENSEMBLE_SCALAR_TYPE
  "Enable use of the ensemble UQ scalar type in stokhos"
  ${Stokhos_ENABLE_Ensemble_Scalar_Type_Default} )

TRIBITS_ADD_OPTION_AND_DEFINE(Stokhos_ENABLE_PCE_Scalar_Type
  HAVE_STOKHOS_PCE_SCALAR_TYPE
  "Enable use of the PCE UQ scalar type in stokhos"
  ${Stokhos_ENABLE_Ensemble_Scalar_Type_Default} )

IF(Stokhos_ENABLE_Ensemble_Scalar_Type AND NOT Stokhos_ENABLE_Sacado)
  MESSAGE(FATAL_ERROR "Ensemble scalar type cannot be enabled unless Sacado is enabled!")
ENDIF()
IF(Stokhos_ENABLE_PCE_Scalar_Type AND NOT Stokhos_ENABLE_Sacado)
  MESSAGE(FATAL_ERROR "PCE scalar type cannot be enabled unless Sacado is enabled!")
ENDIF()

# List of ensemble sizes to explicitly instantiate
IF(Stokhos_ENABLE_Ensemble_Scalar_Type)
  SET(Stokhos_ETI_Ensemble_Sizes "16" CACHE STRING "List of ensemble sizes for explicit instantiation")
  MESSAGE("-- Enabling ensemble sizes ${Stokhos_ETI_Ensemble_Sizes}")
  LIST(LENGTH Stokhos_ETI_Ensemble_Sizes Stokhos_ETI_Num_Ensemble_Sizes)
  IF(Stokhos_ETI_Num_Ensemble_Sizes EQUAL 0)
    MESSAGE(FATAL_ERROR "Must have at least one ensemble size when enabling ensemble scalar type")
  ENDIF()
  # default ensemble size for tests is the first one but user can override
  LIST(GET Stokhos_ETI_Ensemble_Sizes 0 Stokhos_ETI_Default_Ensemble_Size_Default)
  SET(Stokhos_ETI_Default_Ensemble_Size "${Stokhos_ETI_Default_Ensemble_Size_Default}" CACHE STRING "Default ensemble size to use for testing")
  MESSAGE("-- Ensemble size for testing is ${Stokhos_ETI_Default_Ensemble_Size}")
  SET(STOKHOS_DEFAULT_ENSEMBLE_SIZE ${Stokhos_ETI_Default_Ensemble_Size})
  LIST(FIND Stokhos_ETI_Ensemble_Sizes 4 Stokhos_Ensemble_4_Index)
  LIST(FIND Stokhos_ETI_Ensemble_Sizes 8 Stokhos_Ensemble_8_Index)
  LIST(FIND Stokhos_ETI_Ensemble_Sizes 16 Stokhos_Ensemble_16_Index)
  LIST(FIND Stokhos_ETI_Ensemble_Sizes 32 Stokhos_Ensemble_32_Index)
  IF(NOT (Stokhos_Ensemble_4_Index EQUAL -1))
    SET(HAVE_STOKHOS_ENSEMBLE_SIZE_4 1)
  ENDIF()
  IF(NOT (Stokhos_Ensemble_8_Index EQUAL -1))
    SET(HAVE_STOKHOS_ENSEMBLE_SIZE_8 1)
  ENDIF()
  IF(NOT (Stokhos_Ensemble_16_Index EQUAL -1))
    SET(HAVE_STOKHOS_ENSEMBLE_SIZE_16 1)
  ENDIF()
  IF(NOT (Stokhos_Ensemble_32_Index EQUAL -1))
    SET(HAVE_STOKHOS_ENSEMBLE_SIZE_32 1)
  ENDIF()
ENDIF()

# Check if immintrin.h is available
CHECK_INCLUDE_FILE_CXX(immintrin.h HAVE_STOKHOS_INTRINSICS)

#IF(Stokhos_ENABLE_ML)
#  SET(HAVE_STOKHOS_ML ON)
#ENDIF()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()

#
