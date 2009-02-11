
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA AND NOX_ENABLE_LOCA)
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    LOCA_ContinuationManager_XYZT_Example
    SOURCES 
    PeriodicLinearSystem.H
    PeriodicLinearSystem.C
    mainPeriodic.C
    COMM serial mpi
    NUM_MPI_PROCS 2
    ARGS task.xml
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR(
    LOCAContinuationManagerXYZTCopyFiles
    SOURCE_FILES task.xml.in
    DEST_FILES task.xml
    EXEDEPS LOCA_ContinuationManager_XYZT_Example
    )

ENDIF()
