
INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_SOURCE_DIR})
INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR})
INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_SOURCE_DIR}/stk_ngp_test)
INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_SOURCE_DIR}/stk_util)
INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_BINARY_DIR}/stk_util)


IF(Trilinos_ENABLE_STKIO)

FILE(GLOB HEADERS *.h *.hpp)
FILE(GLOB SOURCES *.cpp)

ELSEIF(TPL_ENABLE_MPI)

SET(HEADERS
     ParallelGtestOutput.hpp
     getOption.h
     StkReportRedirector.hpp
     stringAndNumberComparisons.hpp)

SET(SOURCES
     ParallelGtestOutput.cpp
     StkReportRedirector.cpp
     stringAndNumberComparisons.cpp)

ELSE()

SET(HEADERS
    getOption.h
    StkReportRedirector.hpp
    stringAndNumberComparisons.hpp)

SET(SOURCES
    StkReportRedirector.cpp
    stringAndNumberComparisons.cpp)

ENDIF()


IF(Trilinos_ENABLE_STKIO)
  ADD_SUBDIRECTORY(stk_mesh_fixtures)
ENDIF()

INSTALL(FILES ${HEADERS} DESTINATION
    ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_unit_test_utils/)
    
TRIBITS_ADD_LIBRARY(
  stk_unit_test_utils
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  ${Gtest_NO_INSTALL_LIB_OR_HEADERS_ARG}
  )

IF (NOT DEFINED STK_ENABLE_UnitMain)
  SET(STK_ENABLE_UnitMain ON)
ENDIF()

IF (STK_ENABLE_UnitMain)
  FILE(GLOB UNIT_MAIN_SOURCES unit_main_lib/*.cpp)

  TRIBITS_ADD_LIBRARY(
    stk_unit_main
    SOURCES ${UNIT_MAIN_SOURCES}
  )
ENDIF()

