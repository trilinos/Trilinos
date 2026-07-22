INCLUDE(CombinedOption)

#
# A) Define the package
#

TRIBITS_PACKAGE_DECL(Panzer)

#
# B) Set up package-specific options
#

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG
  PANZER_DEBUG
  "Build with debug code enabled."
  ${${PROJECT_NAME}_ENABLE_DEBUG} )

IF (CMAKE_BUILD_TYPE STREQUAL "DEBUG" AND Kokkos_ENABLE_CUDA)
  SET(${PACKAGE_NAME}_ADD_EXPENSIVE_CUDA_TESTS_DEFAULT OFF)
ELSE()
  SET(${PACKAGE_NAME}_ADD_EXPENSIVE_CUDA_TESTS_DEFAULT ON)
ENDIF()
SET(${PACKAGE_NAME}_ADD_EXPENSIVE_CUDA_TESTS
  ${${PACKAGE_NAME}_ADD_EXPENSIVE_CUDA_TESTS_DEFAULT}
  CACHE BOOL
  "Add expensive convergence tests to cuda builds.  Not recommended for CMAKE_BUILD_TYPE=DEBUG and CUDA")
MESSAGE(STATUS "Convergence tests: ${${PACKAGE_NAME}_ADD_EXPENSIVE_CUDA_TESTS}")

SET(${PACKAGE_NAME}_FADTYPE  "Sacado::Fad::DFad<RealType>"
  CACHE STRING
  "Choose the Sacado automatic differentiation scalar type (default is DFad).")

TRIBITS_ADD_ENABLE_TEUCHOS_TIME_MONITOR_OPTION()

GLOBAL_SET(PANZER_UNIT_TEST_MAIN "${PHALANX_UNIT_TEST_MAIN}")


# Explicit Template Instantiation

# ALWAYS FORCE ETI ON FOR PANZER!  ETI in Panzer is for evaluation
# types.  For now we set flags, until we can remove all the ETI code
# for evaluation types.
# TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()
# Copy of above function with hard coded default of ON. The above
# function is hard coded to the Project level default.
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION
  HAVE_${PACKAGE_NAME_UC}_EXPLICIT_INSTANTIATION
  "Enable the use of explicit template instantiation."
  ON
  )
SET(${PARENT_PACKAGE_NAME}_TEMPLATE_DECL_EXT "_decl.hpp")
SET(${PARENT_PACKAGE_NAME}_TEMPLATE_DEF_EXT "_impl.hpp")

TRIBITS_PROCESS_SUBPACKAGES()

TRIBITS_PACKAGE_DEF()

TRIBITS_PACKAGE_POSTPROCESS()
