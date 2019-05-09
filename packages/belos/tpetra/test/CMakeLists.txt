ADD_SUBDIRECTORY(BlockCG)
ADD_SUBDIRECTORY(FixedPoint)
ADD_SUBDIRECTORY(BlockGmres)
ADD_SUBDIRECTORY(LinearSolverFactory)
ADD_SUBDIRECTORY(MultipleSolves)
ADD_SUBDIRECTORY(MVOPTester)
ADD_SUBDIRECTORY(Native)
ADD_SUBDIRECTORY(OrthoManager)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Experimental)
IF(${PACKAGE_NAME}_ENABLE_Experimental)
  # At the moment, all the tests in this subdirectory depend on
  # experimental code.
  ADD_SUBDIRECTORY(BlockGCRODR)
ENDIF() # ${PACKAGE_NAME}_ENABLE_Experimental

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Issue_3235
  SOURCES Issue_3235.cpp
  COMM mpi
  NUM_MPI_PROCS 2
  )
