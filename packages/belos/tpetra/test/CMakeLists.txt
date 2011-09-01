ADD_SUBDIRECTORY(BlockCG)
# ADD_SUBDIRECTORY(BlockGmres)
ADD_SUBDIRECTORY(MVOPTester)
ADD_SUBDIRECTORY(OrthoManager)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Experimental)
IF(${PACKAGE_NAME}_ENABLE_Experimental)
  # 
  # At the moment, all the tests in this subdirectory depend on
  # experimental code.
  # 
  ADD_SUBDIRECTORY(Gmres)
ENDIF() # ${PACKAGE_NAME}_ENABLE_Experimental
