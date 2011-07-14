INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)


INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING
  "${PACKAGE_SOURCE_DIR}/../thyra/core/example/operator_solve")


SET(COMMON_SOURCES test_single_amesos_thyra_solver.cpp)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_amesos_thyra_driver
  SOURCES
    test_amesos_thyra_driver.cpp
    ${COMMON_SOURCES}
  ARGS
    "--show-all-tests"
  COMM serial mpi
  NUM_MPI_PROCS 1
  )


PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyAmesosTestFiles
  DEST_FILES
    bcsstk01.mtx                              
    bcsstk02.mtx                              
    bcsstk04.mtx                              
    nos1.mtx                                      
    nos5.mtx                                      
    FourByFour.mtx                          
    Diagonal.mtx                              
    MissingADiagonal.mtx                              
    KheadK.mtx                                  
    KheadSorted.mtx
    Khead.triS  
  SOURCE_PREFIX "In_"
  SOURCE_DIR "${PACKAGE_SOURCE_DIR}/../amesos/test/Test_Basic"
  EXEDEPS test_amesos_thyra_driver
  )
