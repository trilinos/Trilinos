#
# A) Define the package
#

TRIBITS_PACKAGE(Stokhos)

#
# B) Set up package-specific options
#

TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

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

IF(NOT HAVE_STOKHOS_ENSEMBLE_REDUCT AND NOT Stokhos_ENABLE_CUDA)
  TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Ensemble_GEMV
    HAVE_STOKHOS_ENSEMBLE_GEMV
    "Enable optimized ensemble GEMV implementation."
    ON )
ELSE()
  TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Ensemble_GEMV
    HAVE_STOKHOS_ENSEMBLE_GEMV
    "Enable optimized ensemble GEMV implementation."
    OFF )
ENDIF()

IF(HAVE_STOKHOS_ENSEMBLE_GEMV)
  SET(Stokhos_Ensemble_GEMV_Cache_Size 180224 CACHE STRING "Cache size used in the ensemble GEMV")
  SET(Stokhos_Ensemble_GEMV_Team_Size 1 CACHE STRING "Team size used in the ensemble GEMV")
ENDIF()

SET(Stokhos_ENABLE_Ensemble_Scalar_Type_Default OFF)
SET(Stokhos_ENABLE_PCE_Scalar_Type_Default OFF)
IF(Stokhos_ENABLE_Sacado)
  SET(Stokhos_ENABLE_Ensemble_Scalar_Type_Default ON)
  IF(NOT Stokhos_HIDE_DEPRECATED_CODE)
    SET(Stokhos_ENABLE_PCE_Scalar_Type_Default ON)
  ENDIF()
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(Stokhos_ENABLE_Ensemble_Scalar_Type
  HAVE_STOKHOS_ENSEMBLE_SCALAR_TYPE
  "Enable use of the ensemble UQ scalar type in stokhos"
  ${Stokhos_ENABLE_Ensemble_Scalar_Type_Default} )

TRIBITS_ADD_OPTION_AND_DEFINE(Stokhos_ENABLE_PCE_Scalar_Type
  HAVE_STOKHOS_PCE_SCALAR_TYPE
  "Enable use of the PCE UQ scalar type in stokhos"
  ${Stokhos_ENABLE_PCE_Scalar_Type_Default}  )

IF(Stokhos_ENABLE_PCE_Scalar_Type)
  IF(Stokhos_HIDE_DEPRECATED_CODE)
    MESSAGE(FATAL_ERROR "The PCE scalar type is deprecated in Stokhos but Stokhos_HIDE_DEPRECATED_CODE=ON.  You must set Stokhos_ENABLE_PCE_Scalar_Type=OFF or Stokhos_HIDE_DEPRECATED_CODE=OFF.")
  ELSE()
    MESSAGE(WARNING "The PCE scalar type is deprecated in Stokhos but Stokhos_HIDE_DEPRECATED_CODE=OFF.")
  ENDIF()
ENDIF()

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

  # Make list of enabled ensemble sizes suitable for iteration/searching
  SET(STOKHOS_MPL_ENSEMBLE_SIZES "typedef Sacado::mpl::vector< ")
  FOREACH(ENSEMBLE_SIZE ${Stokhos_ETI_Ensemble_Sizes})
    STRING(APPEND STOKHOS_MPL_ENSEMBLE_SIZES "Sacado::mpl::integral_c<unsigned,${ENSEMBLE_SIZE}>, ")
  ENDFOREACH()
  # Strip last ","
  STRING(LENGTH "${STOKHOS_MPL_ENSEMBLE_SIZES}" ENSEMBLE_LEN)
  MATH(EXPR ENSEMBLE_LEN "${ENSEMBLE_LEN}-2")
  STRING(SUBSTRING "${STOKHOS_MPL_ENSEMBLE_SIZES}" 0 ${ENSEMBLE_LEN} STOKHOS_MPL_ENSEMBLE_SIZES)
  STRING(APPEND STOKHOS_MPL_ENSEMBLE_SIZES " > ETI_Ensemble_Sizes;")
  SET(HAVE_STOKHOS_ENSEMBLE_SIZES "${Stokhos_ETI_Ensemble_Sizes}")
  CONFIGURE_FILE(cmake/Stokhos_Ensemble_Sizes.hpp.in src/Stokhos_Ensemble_Sizes.hpp)
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
