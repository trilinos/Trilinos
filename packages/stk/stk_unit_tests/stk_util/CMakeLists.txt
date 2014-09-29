

ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

if ( TPL_ENABLE_MPI )
  ADD_SUBDIRECTORY(parallel)
ENDIF()

INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_BINARY_DIR}/stk_util)

FILE(GLOB SOURCES diag/*.cpp environment/*.cpp util/*.cpp)
#removing due to dependency on Sierra's unitTestUtils/comparison/stringAndNumberComparisons.h
LIST(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/diag/TimerHowTo.cpp)

APPEND_SET( SOURCES
 UnitTestBoost.cpp
 UtilUnitTestMain.cpp
)

if ( TPL_ENABLE_MPI )
  APPEND_SET( SOURCES
    UnitTestParallel.cpp
  )
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    util_UnitTest
    SOURCES ${SOURCES}
    DEPLIBS stk_util_parallel stk_util_diag
    COMM serial mpi
    NUM_MPI_PROCS 1
    )


