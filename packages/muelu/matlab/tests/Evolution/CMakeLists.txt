TRIBITS_COPY_FILES_TO_BINARY_DIR(Matlab_Evolution_Tests_cp
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  SOURCE_FILES evoEasyDriver.m  evolutionSoCFactory.m laplacianfun.m  matlabParamsEasy.xml classicdrop.xml classicEasyDriver.m Aani.mtx Aiso.mtx Arot.mtx
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  )
