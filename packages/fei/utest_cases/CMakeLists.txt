
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../base/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../support-Trilinos/)

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

#
# Core headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
fei_unit_CommUtils.hpp
fei_unit_CSRMat_CSVec.hpp
fei_unit_DirBC.hpp
fei_unit_impl_utils.hpp
fei_unit_MatrixTraits.hpp
fei_unit_Params.hpp
fei_unit_Reducer.hpp
fei_unit_register_tests.hpp
fei_unit_testcase.hpp
fei_unit_testcontainer.hpp
fei_unit_test_runner.hpp
fei_unit_utils.hpp
  )

APPEND_SET(SOURCES
fei_unit_CommUtils.cpp
fei_unit_CSRMat_CSVec.cpp
fei_unit_DirBC.cpp
fei_unit_impl_utils.cpp
fei_unit_MatrixTraits.cpp
fei_unit_Params.cpp
fei_unit_Reducer.cpp
fei_unit_register_tests.cpp
fei_unit_testcase.cpp
fei_unit_testcontainer.cpp
fei_unit_test_runner.cpp
fei_unit_utils.cpp
  )

IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  fei_utest_cases
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
