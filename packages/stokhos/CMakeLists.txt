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
