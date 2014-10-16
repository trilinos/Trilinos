

IF(Kokkos_ENABLE_Cuda)
  #
  # Add subpackage include path dependencies for CUDA_COMPILE().
  CUDA_INCLUDE_DIRECTORIES( "." )
  CUDA_INCLUDE_DIRECTORIES( ${Tpetra_INCLUDE_DIRS} )
  CUDA_COMPILE(IE_UT_OBJECTS_CUDA ImportExport_UnitTests_CU.cu)
  CUDA_COMPILE(I_U_OBJECTS_CUDA Import_Union_CU.cu)
  CUDA_COMPILE(SMET_OBJECTS_CUDA SubmapExport_CU.cu)
  CUDA_COMPILE(SMIT_OBJECTS_CUDA SubmapImport_CU.cu)
  CUDA_COMPILE(IB5430T_OBJECTS_CUDA ImportBug5430Test_CU.cu)
  CUDA_COMPILE(ETSGCM_OBJECTS_CUDA ExportToStaticGraphCrsMatrix_CU.cu)
  CUDA_COMPILE(MAKEOPTCOLMAP_OBJECTS_CUDA MakeOptColMap_CU.cu)
  CUDA_COMPILE(BUG6170_OBJECTS_CUDA Bug6170_CU.cu)
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ImportExport_UnitTests
  SOURCES
    ImportExport_UnitTests
    ${IE_UT_OBJECTS_CUDA}
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Import_Union
  SOURCES
    Import_Union
    ${I_U_OBJECTS_CUDA}
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SubmapExportTests
  SOURCES
    SubmapExport
    ${SMET_OBJECTS_CUDA}
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SubmapImportTests
  SOURCES
    SubmapImport
    ${SMIT_OBJECTS_CUDA}
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

IF (${PROJECT_NAME}_ENABLE_Epetra)
  IF(NOT Trilinos_NO_32BIT_GLOBAL_INDICES)
  # Tpetra bug 5430:
  # https://software.sandia.gov/bugzilla/show_bug.cgi?id=5430
  #
  # The bug only manifests with 4 MPI processes.  The test currently
  # requires Epetra; it compares Epetra and Tpetra results.
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ImportBug5430
    SOURCES
      ImportBug5430Test.cpp
      ${IB5430T_OBJECTS_CUDA}
    COMM mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 4
    )
  ENDIF ()
ENDIF ()

# Test whether Export to a CrsMatrix with a static graph works.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ExportToStaticGraphCrsMatrix
  SOURCES
    ExportToStaticGraphCrsMatrix.cpp
    ${ETSGCM_OBJECTS_CUDA}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 4
  )

# Test Tpetra::Details::makeOptimizedColMap and
# Tpetra::Details::makeOptimizedColMapAndImport.
#
# This test really only makes sense with an MPI build, since
# optimizing the column Map involves reordering remote indices.  The
# test already exercises the case where there are no remote indices.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MakeOptColMap
  SOURCES
    MakeOptColMap.cpp
    ${MAKEOPTCOLMAP_OBJECTS_CUDA}
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 4
  )

# mfh 18 Jul 2014: Bug 6170 test.  Thanks to Clark Dohrmann for
# contributing this test.  I modified it slightly to avoid file I/O
# and to catch and print exception messages on all processes.  The
# test must be run with exactly 8 MPI processes.
#
# NOTE: The default behavior of the Trilinos test framework is to
# exclude tests that require more than 4 MPI processes.  Thus, this
# test will NOT run by default.  To change this, set the CMake option
# MPI_EXEC_MAX_NUMPROCS to 8 or more in your Trilinos configuration.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Bug6170
  SOURCES
    Bug6170.cpp
    ${BUG6170_OBJECTS_CUDA}
  COMM mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 8
  )
