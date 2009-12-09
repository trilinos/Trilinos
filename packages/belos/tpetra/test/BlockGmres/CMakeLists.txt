
INCLUDE(PackageAddExecutableAndTest)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)

   PACKAGE_ADD_EXECUTABLE_AND_TEST(
     Tpetra_BlockGMRES_hb_DoubleFloat_test
     SOURCES test_bl_gmres_hb_df.cpp 
     ARGS 
     COMM serial mpi
     )
   
   PACKAGE_COPY_FILES_TO_BINARY_DIR(Tpetra_BlockGMRES_hb_DoubleFloat_CopyFiles
     SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/tpetra/test/BlockCG
     SOURCE_FILES bcsstk17.rsa bcsstk18.rsa
     EXEDEPS Tpetra_BlockGMRES_hb_DoubleFloat_test
     )

ENDIF()
