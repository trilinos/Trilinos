
ADD_SUBDIRECTORY(BlockCG)
ADD_SUBDIRECTORY(BlockGmres)

# mfh 10 Jan 2013: MINRES test is taking a long time and causing test timeouts.
# It's not obviously broken, so I'm disabling the test for now until we can find
# a different test problem that doesn't take so long.
#ADD_SUBDIRECTORY(MINRES)
ADD_SUBDIRECTORY(MVOPTester)
IF(${PACKAGE_NAME}_ENABLE_Experimental)
  ADD_SUBDIRECTORY(ProjectedLeastSquaresSolver)
ENDIF()
ADD_SUBDIRECTORY(TFQMR)
