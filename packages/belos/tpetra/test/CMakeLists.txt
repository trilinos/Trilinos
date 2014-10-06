ADD_SUBDIRECTORY(BlockCG)
# ADD_SUBDIRECTORY(BlockGmres)
ADD_SUBDIRECTORY(MVOPTester)
ADD_SUBDIRECTORY(OrthoManager)

ASSERT_DEFINED(Belos_ENABLE_Experimental)
IF(Belos_ENABLE_Experimental)
  # At the moment, all the tests in this subdirectory depend on
  # experimental code.
  ADD_SUBDIRECTORY(Gmres)

  # At the moment, all the tests in this subdirectory depend on
  # experimental code.
  ADD_SUBDIRECTORY(BlockGCRODR)
ENDIF() # ${PACKAGE_NAME}_ENABLE_Experimental

