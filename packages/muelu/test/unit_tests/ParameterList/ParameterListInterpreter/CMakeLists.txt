INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_COPY_FILES_TO_BINARY_DIR(ParameterList_ParameterListInterpreter_cp
  SOURCE_FILES 
  Smoothed-Aggregation.xml   # duplication of factories
  Smoothed-Aggregation2.xml  # use <factories>
  Smoothed-Aggregation3.xml  # test factories as input parameters
  )
