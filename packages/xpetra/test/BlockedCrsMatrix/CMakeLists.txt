INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

#TODO: re-enable
# IF (${PACKAGE_NAME}_ENABLE_Epetra)
#   TRIBITS_ADD_EXECUTABLE_AND_TEST(
#     BlockedCrsMatrix_UnitTests
#     SOURCES
#       BlockedCrsMatrix_UnitTests
#       ${TEUCHOS_STD_UNIT_TEST_MAIN}
#     COMM serial mpi
#     STANDARD_PASS_OUTPUT
#     )
  
#   TRIBITS_COPY_FILES_TO_BINARY_DIR(BlockedCrsMatrix_UnitTests_copyfiles
#     SOURCE_FILES 
#       nsjac_test.mm
#       nslhs_test.mm
#       nsrhs_test.mm
#       A.mat
#       P.mat
#     )
# ENDIF()
