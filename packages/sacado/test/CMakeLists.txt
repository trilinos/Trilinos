
ADD_SUBDIRECTORY(utils)

IF (TPL_ENABLE_CppUnit)
  ADD_SUBDIRECTORY(TestSuite)
ENDIF()

# Not sure how to do this in cmake yet...
##ADD_SUBDIRECTORY(tradoptest)
