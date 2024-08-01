



#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")

# Set the base Trilinos directory in order to find Trilinos_version.h

ASSERT_DEFINED(PROJECT_SOURCE_DIR)

TRIBITS_INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})

# Set headers (for IDEs and installing)

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

# NOTE: Above, we add the XXX_config.h file so that MS Visual Studio
# will know to include this file in the dependency analysis

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../Include)
TRIBITS_INCLUDE_DIRECTORIES(${Trilinos_SOURCE_DIR}/packages/common/auxiliarySoftware/SuiteSparse/src/AMD/Include)
TRIBITS_INCLUDE_DIRECTORIES(${Trilinos_SOURCE_DIR}/packages/common/auxiliarySoftware/SuiteSparse/src/COLAMD/Include)
TRIBITS_INCLUDE_DIRECTORIES(${Trilinos_SOURCE_DIR}/packages/common/auxiliarySoftware/SuiteSparse/src/BTF/Include)
TRIBITS_INCLUDE_DIRECTORIES(${Trilinos_SOURCE_DIR}/packages/common/auxiliarySoftware/SuiteSparse/src/UFconfig)

APPEND_SET(HEADERS
klu2_analyze_given.hpp
klu2_factor.hpp
klu2_refactor.hpp
klu2_analyze.hpp
klu2_free_numeric.hpp 
klu2_scale.hpp
klu2_defaults.hpp
klu2_free_symbolic.hpp
klu2_solve.hpp
klu2_diagnostics.hpp 
klu2.hpp
klu2_sort.hpp
klu2_dump.hpp
klu2_kernel.hpp
klu2_tsolve.hpp
klu2_extract.hpp
klu2_memory.hpp
  )

#
# C) Define the targets for package's test
#

TRIBITS_ADD_TEST_DIRECTORIES(tests)
