TRIBITS_SUBPACKAGE(Unit_test_utils)

INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_SOURCE_DIR})

FILE(GLOB HEADERS *.hpp)
#disabling due to creation of circular dependency on stk_unit_tests mesh_fixture
LIST(REMOVE_ITEM HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/unittestMeshUtils.hpp")
FILE(GLOB SOURCES *.cpp)
#disabling due to creation of circular dependency on stk_unit_tests mesh_fixture
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/unittestMeshUtils.cpp")

TRIBITS_ADD_LIBRARY(
  stk_unit_test_utils
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  #NO_INSTALL_LIB_OR_HEADERS
  #TESTONLY
  )

TRIBITS_SUBPACKAGE_POSTPROCESS()
