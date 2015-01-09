

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)

   TRIBITS_ADD_EXECUTABLE_AND_TEST(
     Tpetra_BlockGMRES_hb_DoubleFloat_test
     SOURCES test_bl_gmres_hb_df.cpp 
     ARGS 
     COMM serial mpi
     )
   
   TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_BlockGMRES_hb_DoubleFloat_CopyFiles
     SOURCE_DIR ${Belos_SOURCE_DIR}/tpetra/test/BlockCG
     SOURCE_FILES bcsstk17.rsa bcsstk18.rsa
     EXEDEPS Tpetra_BlockGMRES_hb_DoubleFloat_test
     )

ENDIF()
