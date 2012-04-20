INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${TEUCHOS_STD_UNIT_TEST_MAIN})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tEpetraLinearObjFactory
  SOURCES tEpetraLinearObjFactory.cpp UnitTest_UniqueGlobalIndexer.cpp ${UNIT_TEST_DRIVER}
  COMM mpi
  NUM_MPI_PROCS 2
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tTpetraLinearObjFactory
  SOURCES tTpetraLinearObjFactory.cpp UnitTest_UniqueGlobalIndexer.cpp ${UNIT_TEST_DRIVER}
  COMM mpi
  NUM_MPI_PROCS 2
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tBlockedLinearObjFactory
  SOURCES tBlockedLinearObjFactory.cpp UnitTest_UniqueGlobalIndexer.cpp ${UNIT_TEST_DRIVER}
  COMM mpi
  NUM_MPI_PROCS 2
  )

IF(${PACKAGE_NAME}_ENABLE_Stokhos)
   TRIBITS_ADD_EXECUTABLE_AND_TEST(
     tSGEpetraLinearObjFactory
     SOURCES tSGEpetraLinearObjFactory.cpp UnitTest_UniqueGlobalIndexer.cpp ${UNIT_TEST_DRIVER}
     COMM mpi
     NUM_MPI_PROCS 2
     )
ENDIF()
