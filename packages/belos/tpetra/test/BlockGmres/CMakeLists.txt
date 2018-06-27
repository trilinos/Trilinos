

ASSERT_DEFINED(Anasazi_SOURCE_DIR)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)

IF (${PACKAGE_NAME}_ENABLE_Triutils)

   TRIBITS_ADD_EXECUTABLE_AND_TEST(
     Tpetra_BlockGMRES_hb_DoubleFloat_test
     SOURCES test_bl_gmres_hb_df.cpp 
     ARGS "--verbose"
     COMM serial mpi
     )
   
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_BlockGMRES_hb_test
  SOURCES test_bl_gmres_hb.cpp
  ARGS "--verbose"
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_BlockGMRES_hb_CopyFiles
  SOURCE_DIR ${Anasazi_SOURCE_DIR}/testmatrices
  SOURCE_FILES bcsstk14.hb
  EXEDEPS Tpetra_BlockGMRES_hb_test
  )

