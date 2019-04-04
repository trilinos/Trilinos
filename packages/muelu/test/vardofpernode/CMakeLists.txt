INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters/belos)

IF (${PACKAGE_NAME}_ENABLE_Belos)

TRIBITS_ADD_EXECUTABLE(
    VarDofDriver
    SOURCES VarDofDriver.cpp
    COMM serial mpi
  )

ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(VarDofDriver_cp
  SOURCE_FILES xxx.mm yyy.mm driver.xml theMatrix.m rhs.mm 1ProcLinear0 1proc0 1dofPresent0 2ProcLinear0  2ProcLinear1  2dofPresent0  2dofPresent1  2proc0  2proc1 7ProcLinear0 7ProcLinear1 7ProcLinear2 7ProcLinear3 7ProcLinear4 7ProcLinear5 7ProcLinear6 7dofPresent0 7dofPresent1 7dofPresent2 7dofPresent3 7dofPresent4 7dofPresent5 7dofPresent6 7proc0 7proc1 7proc2 7proc3 7proc4 7proc5 7proc6

)

IF (${PACKAGE_NAME}_ENABLE_Belos)
      TRIBITS_ADD_TEST(
        VarDofDriver
        NAME "VarDofDriver"
        ARGS "--xml=driver.xml"
        NUM_MPI_PROCS 1
        COMM serial mpi
        )

      TRIBITS_ADD_TEST(
        VarDofDriver
        NAME "VarDofDriver"
        ARGS "--xml=driver.xml"
        NUM_MPI_PROCS 2
        COMM mpi # HAVE_MPI required
        )

      TRIBITS_ADD_TEST(
        VarDofDriver
        NAME "VarDofDriver"
        ARGS "--xml=driver.xml"
        NUM_MPI_PROCS 7
        COMM mpi # HAVE_MPI required
        )
ENDIF()


