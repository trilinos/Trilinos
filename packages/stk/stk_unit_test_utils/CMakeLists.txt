TRIBITS_SUBPACKAGE(Unit_test_utils)

ADD_SUBDIRECTORY(comparison)

INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_SOURCE_DIR})

FILE(GLOB HEADERS *.hpp)
#disabling due to creation of circular dependency on stk_unit_tests mesh_fixture
LIST(REMOVE_ITEM HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/unittestMeshUtils.hpp")
FILE(GLOB SOURCES *.cpp)
#disabling due to creation of circular dependency on stk_unit_tests mesh_fixture
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/unittestMeshUtils.cpp")

INSTALL(FILES ${HEADERS} DESTINATION
    ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_unit_test_utils/)
    
TRIBITS_ADD_LIBRARY(
  stk_unit_test_utils
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  #NO_INSTALL_LIB_OR_HEADERS
  #TESTONLY
  )

IF (NOT DEFINED STK_ENABLE_UnitMain)
  SET(STK_ENABLE_UnitMain ON)
ENDIF()

IF (STK_ENABLE_UnitMain)
  FILE(GLOB UNIT_MAIN_SOURCES unit_main_lib/*.cpp)

  TRIBITS_ADD_LIBRARY(
    stk_unit_main
    SOURCES ${UNIT_MAIN_SOURCES}
    DEPLIBS stk_unit_test_utils
  )
ENDIF()

TRIBITS_SUBPACKAGE_POSTPROCESS()
