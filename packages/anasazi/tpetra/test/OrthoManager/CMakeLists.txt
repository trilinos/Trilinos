
INCLUDE(PackageAddExecutableAndTest)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)

  # This test depends on complex arithmetic support.
  IF (${PACKAGE_NAME}_ENABLE_COMPLEX)
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      Tpetra_OrthoManager_Complex_test
      SOURCES cxx_main_complex.cpp 
      ARGS 
      COMM serial mpi
      )

    PACKAGE_COPY_FILES_TO_BINARY_DIR(Tpetra_OrthoManager_ComplexCopyFiles
      SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/anasazi/testmatrices
      SOURCE_FILES mhd1280b.cua
      EXEDEPS Tpetra_OrthoManager_Complex_test
      )
  ENDIF ()

  IF (${PACKAGE_NAME}_ENABLE_TSQR)
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      Tpetra_OrthoManager_Real_Tsqr_test
      SOURCES cxx_main_tsqr.cpp 
      ARGS 
      COMM serial mpi
      )
  ENDIF ()
ENDIF()
