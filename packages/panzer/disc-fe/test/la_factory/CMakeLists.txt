
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

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
INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING
                    ${CMAKE_CURRENT_SOURCE_DIR}
                    ${CMAKE_CURRENT_SOURCE_DIR}/../../../dof-mgr/test/dofmngr_test)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  LinearObjFactory_Tests
  SOURCES tTpetraLinearObjFactory.cpp 
          tBlockedLinearObjFactory.cpp 
          tBlockedTpetraLinearObjFactory.cpp 
          tTpetra_GlbEvalData.cpp 
          tEpetra_GlbEvalData.cpp 
          UnitTest_GlobalIndexer.cpp 
          ${CMAKE_CURRENT_SOURCE_DIR}/../../../dof-mgr/test/dofmngr_test/UnitTest_ConnManager.cpp
          ${UNIT_TEST_DRIVER}
  COMM mpi
  NUM_MPI_PROCS 2
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tCloneLOF
  SOURCES tCloneLOF.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../../../dof-mgr/test/dofmngr_test/UnitTest_ConnManager.cpp ${UNIT_TEST_DRIVER}
  COMM mpi
  NUM_MPI_PROCS 2
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tEpetra_LOF_FilteredGlobalIndexer
  SOURCES tEpetra_LOF_FilteredGlobalIndexer.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../../../dof-mgr/test/dofmngr_test/UnitTest_ConnManager.cpp ${UNIT_TEST_DRIVER}
  COMM mpi
  NUM_MPI_PROCS 2
  )

