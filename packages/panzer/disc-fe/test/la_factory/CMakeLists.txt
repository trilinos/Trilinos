
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${PANZER_UNIT_TEST_MAIN})

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   tTpetraLinearObjFactory
#   SOURCES tTpetraLinearObjFactory.cpp UnitTest_GlobalIndexer.cpp ${UNIT_TEST_DRIVER}
#   COMM mpi
#   NUM_MPI_PROCS 2
#   )
#
# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   tBlockedLinearObjFactory
#   SOURCES tBlockedLinearObjFactory.cpp UnitTest_GlobalIndexer.cpp
#           ${CMAKE_CURRENT_SOURCE_DIR}/../dofmngr_test/UnitTest_ConnManager.cpp
#           ${UNIT_TEST_DRIVER}
#   COMM mpi
#   NUM_MPI_PROCS 2
#   )
#
# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   tBlockedTpetraLinearObjFactory
#   SOURCES tBlockedTpetraLinearObjFactory.cpp UnitTest_GlobalIndexer.cpp ${UNIT_TEST_DRIVER}
#   COMM mpi
#   NUM_MPI_PROCS 2
#   )
#
# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   tEpetra_GlbEvalData
#   SOURCES tEpetra_GlbEvalData.cpp UnitTest_GlobalIndexer.cpp ${UNIT_TEST_DRIVER}
#   COMM mpi
#   NUM_MPI_PROCS 2
#   )

# The intention is that this INCLUDE only happens for tests below this point.
# If this assumption is incorrect, or the build system is changed this might
# have to be fixed!
###############################################################################
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING
                    ${CMAKE_CURRENT_SOURCE_DIR}
                    ${CMAKE_CURRENT_SOURCE_DIR}/../../../dof-mgr/test/dofmngr_test)

SET(FILES tTpetraLinearObjFactory.cpp tTpetra_GlbEvalData.cpp UnitTest_GlobalIndexer.cpp tBlockedTpetraLinearObjFactory.cpp tTpetraFEGraphAssembly.cpp)

IF(PANZER_HAVE_EPETRA_STACK)
  APPEND_SET(FILES tBlockedLinearObjFactory.cpp
             tEpetra_GlbEvalData.cpp)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    tCloneLOF_Epetra
    SOURCES tEpetra_tCloneLOF.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../../../dof-mgr/test/dofmngr_test/UnitTest_ConnManager.cpp ${UNIT_TEST_DRIVER}
    COMM mpi
    NUM_MPI_PROCS 2
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    LOF_FilteredGlobalIndexer_Epetra
    SOURCES tEpetra_LOF_FilteredGlobalIndexer.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../../../dof-mgr/test/dofmngr_test/UnitTest_ConnManager.cpp ${UNIT_TEST_DRIVER}
    COMM mpi
    NUM_MPI_PROCS 2
  )
ENDIF(PANZER_HAVE_EPETRA_STACK)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  LinearObjFactory_Tests
  SOURCES ${FILES}
          ${CMAKE_CURRENT_SOURCE_DIR}/../../../dof-mgr/test/dofmngr_test/UnitTest_ConnManager.cpp
          ${UNIT_TEST_DRIVER}
  COMM mpi
  NUM_MPI_PROCS 2
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tCloneLOF_Tpetra
  SOURCES tTpetra_tCloneLOF.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../../../dof-mgr/test/dofmngr_test/UnitTest_ConnManager.cpp ${UNIT_TEST_DRIVER}
  COMM mpi
  NUM_MPI_PROCS 2
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  LOF_FilteredGlobalIndexer_Tpetra
  SOURCES tTpetra_LOF_FilteredGlobalIndexer.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../../../dof-mgr/test/dofmngr_test/UnitTest_ConnManager.cpp ${UNIT_TEST_DRIVER}
  COMM mpi
  NUM_MPI_PROCS 2
)

# Time (and compare) building the FE graph on device vs. by insertion, on a
# structured hex mesh.  The test runs a small mesh and checks the graphs are identical.
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../dof-mgr/test/cartesian_topology)
TRIBITS_ADD_EXECUTABLE(
  FEGraphAssemblyBenchmark
  SOURCES FEGraphAssemblyBenchmark.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/../../../dof-mgr/test/cartesian_topology/CartesianConnManager.cpp
  COMM mpi
  )

TRIBITS_ADD_TEST(
  FEGraphAssemblyBenchmark
  NAME FEGraphAssembly_Identical
  ARGS "--nx=6 --ny=5 --nz=4 --px=2 --py=2 --pz=1 --repetitions=1"
  COMM mpi
  NUM_MPI_PROCS 4
  PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
  )
