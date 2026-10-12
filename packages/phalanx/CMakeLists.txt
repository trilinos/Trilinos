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

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_SHARED_SPACE
  PHX_ENABLE_SHARED_SPACE
  "If set to ON, the default memory space for phalanx will be Kokkos::SharedSpace (default=OFF)"
  OFF )
SET(${PACKAGE_NAME}_DEFAULT_EXECUTION_SPACE "Kokkos::DefaultExecutionSpace" CACHE STRING "The default execution space")

# The memory space is derived, never chosen: either the execution space's own,
# or Kokkos::SharedSpace.
IF (${PACKAGE_NAME}_ENABLE_SHARED_SPACE)
  SET(PHX_DEFAULT_MEMORY_SPACE_CXX "Kokkos::SharedSpace")
ELSE()
  SET(PHX_DEFAULT_MEMORY_SPACE_CXX "${${PACKAGE_NAME}_DEFAULT_EXECUTION_SPACE}::memory_space")
ENDIF()

# Whether PHX::Device's memory space is a DIFFERENT TYPE from its execution
# space's own.  This decides whether PHX::Device is already among the
# per-backend Intrepid2 instantiations in Panzer_IntrepidBasisFactory.cpp, so
# it has to be exact in both directions: wrong one way gives undefined
# references at every call site, wrong the other gives a duplicate explicit
# instantiation.  Both cases below are decided by Kokkos' own rules rather than
# by comparing type names, which could not tell one type spelled two ways from
# two types.
#
#   shared space off   The memory space IS the execution space's own.
#   shared space on    Kokkos::SharedSpace is a distinct space only when a GPU
#                      backend provides it; with host backends alone Kokkos
#                      defines it AS Kokkos::HostSpace, which is the host
#                      execution space's own.
IF (${PACKAGE_NAME}_ENABLE_SHARED_SPACE AND
    (Kokkos_ENABLE_CUDA OR Kokkos_ENABLE_HIP OR Kokkos_ENABLE_SYCL))
  SET(PHX_MEMORY_SPACE_IS_NOT_EXECUTION_SPACE_MEMORY_SPACE ON)
ELSE()
  SET(PHX_MEMORY_SPACE_IS_NOT_EXECUTION_SPACE_MEMORY_SPACE OFF)
ENDIF()

SET(${PACKAGE_NAME}_INDEX_SIZE_TYPE "int" CACHE STRING "The default index type, choose a C++ integer type or PHX::DefaultExecutionSpace::size_type (default=int)")

# PHX::Device is a Kokkos::Device, which carries an execution space and a
# memory space independently.  It used to be a bare execution space, so code
# that still writes e.g. Kokkos::RangePolicy<PHX::Device> needs this option
# until it migrates.
TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEPRECATED_DEVICE_AS_EXECUTION_SPACE
  PHX_DEPRECATED_DEVICE_AS_EXECUTION_SPACE
  "DEPRECATED. Keep PHX::Device as an execution space rather than a Kokkos::Device, for applications that have not migrated yet. Shared space support requires this to be OFF."
  OFF )
IF (${PACKAGE_NAME}_ENABLE_DEPRECATED_DEVICE_AS_EXECUTION_SPACE)
  MESSAGE(DEPRECATION "Phalanx: PHX::Device as an execution space is deprecated. Migrate uses of PHX::Device in execution space positions to PHX::ExecutionSpace -- packages/phalanx/scripts/migrate_phx_device.py does most of it -- then reconfigure with ${PACKAGE_NAME}_ENABLE_DEPRECATED_DEVICE_AS_EXECUTION_SPACE=OFF.")
ENDIF()

MESSAGE(STATUS "Enable Shared Space: ${${PACKAGE_NAME}_ENABLE_SHARED_SPACE}")
MESSAGE(STATUS "Default Execution Space: ${${PACKAGE_NAME}_DEFAULT_EXECUTION_SPACE}")
MESSAGE(STATUS "Default Memory Space: ${PHX_DEFAULT_MEMORY_SPACE_CXX}")

# Map the keywords this option used to take onto the C++ types they selected,
# so existing configure scripts and caches keep working.
SET(PHX_LEGACY_INDEX_SIZE_TYPES KOKKOS INT UINT LONGINT ULONGINT)
# PHX::DefaultExecutionSpace, not PHX::ExecutionSpace: this type is
# substituted into Phalanx_config.hpp, which declares DefaultExecutionSpace
# itself but is included BY the header that declares ExecutionSpace.
SET(PHX_LEGACY_INDEX_SIZE_CXX_TYPES
  "PHX::DefaultExecutionSpace::size_type" "int" "unsigned int" "long int" "unsigned long int")
LIST(FIND PHX_LEGACY_INDEX_SIZE_TYPES "${${PACKAGE_NAME}_INDEX_SIZE_TYPE}" PHX_LEGACY_INDEX)
IF (PHX_LEGACY_INDEX GREATER -1)
  LIST(GET PHX_LEGACY_INDEX_SIZE_CXX_TYPES ${PHX_LEGACY_INDEX} PHX_INDEX_SIZE_TYPE_CXX)
  MESSAGE(DEPRECATION "${PACKAGE_NAME}_INDEX_SIZE_TYPE=${${PACKAGE_NAME}_INDEX_SIZE_TYPE} is deprecated; set it to the C++ type instead: \"${PHX_INDEX_SIZE_TYPE_CXX}\"")
ELSE()
  SET(PHX_INDEX_SIZE_TYPE_CXX "${${PACKAGE_NAME}_INDEX_SIZE_TYPE}")
ENDIF()
MESSAGE(STATUS "Index Size Type: ${PHX_INDEX_SIZE_TYPE_CXX}")

IF (${PACKAGE_NAME}_ENABLE_SHARED_SPACE AND
    ${PACKAGE_NAME}_ENABLE_DEPRECATED_DEVICE_AS_EXECUTION_SPACE)
  MESSAGE(FATAL_ERROR "ERROR (Phalanx): ${PACKAGE_NAME}_ENABLE_SHARED_SPACE needs PHX::Device to be a Kokkos::Device so that the configured memory space is honored. It cannot be used with ${PACKAGE_NAME}_ENABLE_DEPRECATED_DEVICE_AS_EXECUTION_SPACE=ON, which makes PHX::Device an execution space again and so forces the memory space to that execution space's own.")
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

# Calling a virtual function from device code is not supported on SYCL
# without the experimental sycl_ext_oneapi_virtual_functions extension, so the
# tests that do it are not built on a SYCL build unless this is turned on.
# Other backends are unaffected.
OPTION(${PACKAGE_NAME}_ENABLE_VIRTUAL_FUNCTIONS_ON_SYCL_DEVICE
  "Build the tests that call virtual functions from device code on a SYCL build.  Requires the experimental sycl_ext_oneapi_virtual_functions extension; Phalanx adds -fsycl-allow-virtual-functions to those tests (default=OFF)."
  OFF )
IF (Kokkos_ENABLE_SYCL)
  MESSAGE(STATUS "Virtual functions on SYCL device: ${${PACKAGE_NAME}_ENABLE_VIRTUAL_FUNCTIONS_ON_SYCL_DEVICE}")
ENDIF()

# Experimental device dag support
TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEVICE_DAG
  PHX_ENABLE_DEVICE_DAG
  "The DAG is traversed on device from a single parallel_for call. Requires relocatable device code for CUDA builds."
  OFF )
MESSAGE(STATUS "Experimental Device DAG support: ${PHX_ENABLE_DEVICE_DAG}")

# Enforce that if device type is cuda and device dag enabled then relocatable device code must be enabled
IF(${PACKAGE_NAME}_ENABLE_DEVICE_DAG)
  IF(Kokkos_ENABLE_CUDA)
    IF(NOT Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)
      MESSAGE(FATAL_ERROR " Enabling ${PACKAGE_NAME}_ENABLE_DEVICE_DAG support with a CUDA device type in Phalanx requires compilation with relocatable device code! Please reconfigure with Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=ON.")
    ENDIF()
  ENDIF()
  # The device DAG dispatches through PHX::DeviceEvaluator's virtual functions
  # inside a kernel, which SYCL only supports through the experimental
  # sycl_ext_oneapi_virtual_functions extension.
  IF(Kokkos_ENABLE_SYCL)
    IF(NOT ${PACKAGE_NAME}_ENABLE_VIRTUAL_FUNCTIONS_ON_SYCL_DEVICE)
      MESSAGE(FATAL_ERROR "Phalanx configure ERROR: ${PACKAGE_NAME}_ENABLE_DEVICE_DAG calls virtual functions from device code, which SYCL does not support. Reconfigure with ${PACKAGE_NAME}_ENABLE_DEVICE_DAG=OFF, or set ${PACKAGE_NAME}_ENABLE_VIRTUAL_FUNCTIONS_ON_SYCL_DEVICE=ON to try it with the experimental sycl_ext_oneapi_virtual_functions extension.")
    ELSE()
      MESSAGE(WARNING "Phalanx: ${PACKAGE_NAME}_ENABLE_DEVICE_DAG on SYCL needs -fsycl-allow-virtual-functions on every translation unit that compiles a DeviceEvaluator, including in downstream packages, so Phalanx cannot add it for you. Put it in CMAKE_CXX_FLAGS yourself.")
    ENDIF()
  ENDIF()
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_KOKKOS_AMT
  PHX_ENABLE_KOKKOS_AMT
  "Enable the Kokkos Asynchronous Multi-tasking capability. Default is off."
  OFF )
MESSAGE(STATUS "Experimental Kokkos AMT support: ${PHX_ENABLE_KOKKOS_AMT}")
IF(${PACKAGE_NAME}_ENABLE_KOKKOS_AMT)
  # Kokkos picks Threads as the default execution space only when no GPU
  # backend and no OpenMP is enabled, so mirror that order here.
  IF(Kokkos_ENABLE_THREADS AND NOT Kokkos_ENABLE_OPENMP AND
     NOT Kokkos_ENABLE_CUDA AND NOT Kokkos_ENABLE_HIP AND NOT Kokkos_ENABLE_SYCL)
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
