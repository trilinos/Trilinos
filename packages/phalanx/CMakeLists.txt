INCLUDE(CombinedOption)

#
# A) Define the package
#

TRIBITS_PACKAGE(Phalanx)

#
# B) Set up package-specific options
#

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ALLOW_MULTIPLE_EVALUATORS_FOR_SAME_FIELD
  PHX_ALLOW_MULTIPLE_EVALUATORS_FOR_SAME_FIELD
  "Enable to use backwards compatible mode where multiple evaluators can be registered that evalaute the same field. Phalanx will use the first evaluator it finds. Default to OFF but unless Panzer is enabled. For backwards compatibility with Panzer users, we default this to ON."
  "${${PROJECT_NAME}_ENABLE_Panzer}" )
MESSAGE(STATUS "Allow Multiple Evaluator Registration for Same Field: ${${PACKAGE_NAME}_ALLOW_MULTIPLE_EVALUATORS_FOR_SAME_FIELD}")

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_EXPLICIT_TEMPLATE_INSTANTIATION
  PHX_ETI
  "Build with explicit template instation enabled in the examples."
  OFF )

OPTION(${PACKAGE_NAME}_ENABLE_SHARED_SPACE "If set to ON, the default memory space for phalanx will be Kokkos::SharedSpace (default=OFF)" OFF)
SET(${PACKAGE_NAME}_DEFAULT_EXECUTION_SPACE "Kokkos::DefaultExecutionSpace" CACHE STRING "The default execution space")
IF (${PACKAGE_NAME}_ENABLE_SHARED_SPACE)
  SET(${PACKAGE_NAME}_DEFAULT_MEMORY_SPACE "Kokkos::SharedSpace" CACHE STRING "The default memory space")
ELSE()
  SET(${PACKAGE_NAME}_DEFAULT_MEMORY_SPACE "Kokkos::DefaultExecutionSpace::memory_space" CACHE STRING "The default memory space")
ENDIF()
SET(${PACKAGE_NAME}_INDEX_SIZE_TYPE "int" CACHE STRING "The default index type, choose a C++ integer type or PHX::DefaultExecutionSpace::size_type (default=int)")

MESSAGE(STATUS "Enable Shared Space: ${${PACKAGE_NAME}_ENABLE_SHARED_SPACE}")
MESSAGE(STATUS "Default Execution Space: ${${PACKAGE_NAME}_DEFAULT_EXECUTION_SPACE}")
MESSAGE(STATUS "Default Memory Space: ${${PACKAGE_NAME}_DEFAULT_MEMORY_SPACE}")
MESSAGE(STATUS "Index Size Type: ${${PACKAGE_NAME}_INDEX_SIZE_TYPE}")

IF (${${PACKAGE_NAME}_ENABLE_SHARED_SPACE} AND NOT "${${PACKAGE_NAME}_DEFAULT_MEMORY_SPACE}" STREQUAL "Kokkos::SharedSpace")
  MESSAGE(FATAL_ERROR "ERROR (Phalanx): You cannot enable SharedSpace and set the default memory space. Do one or the other.")
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG
  PHX_DEBUG
  "Build with debug code enabled. Typically, this enables array bounds checking and runtime checking of rank accessors in multidimensional arrays."
  OFF )
MESSAGE(STATUS "Enable DEBUG: ${${PACKAGE_NAME}_ENABLE_DEBUG}")

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_TEUCHOS_TIME_MONITOR
  PHX_TEUCHOS_TIME_MONITOR
  "Build with Teuchos TimeMonitors wrapped around each evaluator."
  "${${PROJECT_NAME}_ENABLE_TEUCHOS_TIME_MONITOR}" )

# Use only for CUDA builds
MESSAGE(STATUS "Sacado Hierarchic DFAD support: ${Sacado_ENABLE_HIERARCHICAL_DFAD}")

# Experimental device dag support
TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEVICE_DAG
  PHX_ENABLE_DEVICE_DAG
  "The DAG is traversed on device from a single parallel_for call. Requires relocatable device code for CUDA builds."
  OFF )
MESSAGE(STATUS "Experimental Device DAG support: ${PHX_ENABLE_DEVICE_DAG}")

# Enforce that if device type is cuda and device dag enabled then relocatable device code must be enabled
IF(${PACKAGE_NAME}_ENABLE_DEVICE_DAG)
  IF(PHX_KOKKOS_DEVICE_TYPE_CUDA)
    IF(NOT Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)
      MESSAGE(FATAL_ERROR, " Enabling ${PACKAGE_NAME}_ENABLE_DEVICE_DAG support with a CUDA device type in Phalanx requires compilation with relocatable device code! Please reconfigure with Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=ON.")
    ENDIF()
  ENDIF()
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_KOKKOS_AMT
  PHX_ENABLE_KOKKOS_AMT
  "Enable the Kokkos Asynchronous Multi-tasking capability. Default is off."
  OFF )
MESSAGE(STATUS "Experimental Kokkos AMT support: ${PHX_ENABLE_KOKKOS_AMT}")
IF(${PACKAGE_NAME}_ENABLE_KOKKOS_AMT)
  IF(PHX_KOKKOS_DEVICE_TYPE_THREAD)
  MESSAGE(FATAL_ERROR "Phalanx configure ERROR: Kokkos AMT support is enabled but is not implemented for Pthreads execution space. Please change your execution space or disable AMT support in Phalanx.")
ENDIF()
ENDIF(${PACKAGE_NAME}_ENABLE_KOKKOS_AMT)

# For use by Phalanx and Panzer unit tests
GLOBAL_SET(PHALANX_UNIT_TEST_MAIN
  "${PACKAGE_SOURCE_DIR}/test/Utilities/Phalanx_UnitTestMain.cpp")

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

INSTALL(PROGRAMS scripts/phalanx_create_evaluator.py 
	DESTINATION bin) 

TRIBITS_PACKAGE_POSTPROCESS()
