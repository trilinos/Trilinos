
IF(Kokkos_ENABLE_Cuda)
  #
  # Add subpackage include path dependencies for CUDA_COMPILE().
  CUDA_INCLUDE_DIRECTORIES( "." )
  CUDA_INCLUDE_DIRECTORIES( ${Tpetra_INCLUDE_DIRS} )
  CUDA_COMPILE(M_OBJECTS_CUDA main_CU.cu)
  CUDA_COMPILE(RMT_UT_OBJECTS_CUDA RowMatrixTransposer_UnitTests_CU.cu)
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
	RowMatrixTransposer_test
	SOURCES 
          main.cpp
          ${M_OBJECTS_CUDA}
        )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  RowMatrixTransposer_UnitTests
  SOURCES
    RowMatrixTransposer_UnitTests
    ${RMT_UT_OBJECTS_CUDA}
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(TransposeFileCopies
  SOURCE_FILES a.mtx atrans.mtx
  EXEDEPS RowMatrixTransposer_UnitTests
  )

