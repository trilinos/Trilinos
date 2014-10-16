

SET(ARGS "--filedir=${CMAKE_CURRENT_BINARY_DIR}/")

IF (SITE STREQUAL "gabriel.sandia.gov")
  SET(ARGS
    "${ARGS} --not-unit-test=EpetraRowMatrix_int_ComplexFloat_FullMatrixComplex_UnitTest")
ENDIF()

IF(Kokkos_ENABLE_Cuda)
  #
  # Add subpackage include path dependencies for CUDA_COMPILE().
  CUDA_INCLUDE_DIRECTORIES( "." )
  CUDA_INCLUDE_DIRECTORIES( ${Tpetra_INCLUDE_DIRS} )
  CUDA_COMPILE(EPETRA_RM_UT_OBJECTS_CUDA EpetraRowMatrix_UnitTests_CU.cu)
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  EpetraRowMatrix_UnitTests
  SOURCES
    EpetraRowMatrix_UnitTests
    ${EPETRA_RM_UT_OBJECTS_CUDA}
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  STANDARD_PASS_OUTPUT
  )


TRIBITS_COPY_FILES_TO_BINARY_DIR(TpetraEpetraRowMatrixCopyFiles
  SOURCE_FILES west0067.rua mhd1280b.cua
  EXEDEPS EpetraRowMatrix_UnitTests
  )
