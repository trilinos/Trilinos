INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)


INCLUDE_DIRECTORIES("${PACKAGE_SOURCE_DIR}/../thyra/example/operator_solve")


SET(COMMON_SOURCES test_single_aztecoo_thyra_solver.cpp)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_single_aztecoo_thyra_solver_driver
  SOURCES
    test_single_aztecoo_thyra_solver_driver.cpp
    ${COMMON_SOURCES}
  ARGS
    "--matrix-file=FourByFour.mtx --show-all-tests"
  COMM serial mpi
  NUM_MPI_PROCS 1
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_aztecoo_thyra_driver
  SOURCES
    test_aztecoo_thyra_driver.cpp
    ${COMMON_SOURCES}
  ARGS
    "--show-all-tests"
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

ASSERT_DEFINED(CMAKE_HOST_SYSTEM_NAME)
IF (CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux")
  SET(testEpetraOperatorWrapper_tol 1e-4)
ELSE()
  SET(testEpetraOperatorWrapper_tol 1e-3)
ENDIF()

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  testEpetraOperatorWrapper
  SOURCES
    testEpetraOperatorWrapper.cpp
    ${COMMON_SOURCES}
  ARGS --tol=${testEpetraOperatorWrapper_tol}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )


PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyAztecOOTestFiles
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
  SOURCE_PREFIX
    "In_"
  SOURCE_DIR
    "${PACKAGE_SOURCE_DIR}/../amesos/test/Test_Basic"
  EXEDEPS
    test_single_aztecoo_thyra_solver_driver
    test_aztecoo_thyra_driver
    testEpetraOperatorWrapper
  )
