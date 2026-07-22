TRIBITS_COPY_FILES_TO_BINARY_DIR(
  Matlab_Evolution_Tests_cp
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  SOURCE_FILES
    Aani.mtx
    Aiso.mtx
    Arot.mtx
    classicdrop.xml
    classicEasyDriver.m
    evoEasyDriver.m
    evolutionSoCFactory.m
    laplacianfun.m
    matlabParamsEasy.xml
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)
