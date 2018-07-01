ADD_SUBDIRECTORY(BlockCG)
ADD_SUBDIRECTORY(FixedPoint)
ADD_SUBDIRECTORY(BlockGmres)
ADD_SUBDIRECTORY(LinearSolverFactory)
ADD_SUBDIRECTORY(MultipleSolves)
ADD_SUBDIRECTORY(MVOPTester)
ADD_SUBDIRECTORY(OrthoManager)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Experimental)
IF(${PACKAGE_NAME}_ENABLE_Experimental)
  # At the moment, all the tests in this subdirectory depend on
  # experimental code.
  ADD_SUBDIRECTORY(BlockGCRODR)
ENDIF() # ${PACKAGE_NAME}_ENABLE_Experimental

