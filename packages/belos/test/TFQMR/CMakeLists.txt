
INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddTest)
INCLUDE(PackageCopyFilesToBinaryDir)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
    
  PACKAGE_ADD_EXECUTABLE(
    tfqmr_hb
    SOURCES test_tfqmr_hb.cpp createEpetraProblem.cpp
    COMM serial mpi
  )
  
  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
  IF (${PACKAGE_NAME}_ENABLE_Ifpack)

  PACKAGE_ADD_EXECUTABLE(
    ptfqmr_hb
    SOURCES test_ptfqmr_hb.cpp createEpetraProblem.cpp
    COMM serial mpi
  )

  ENDIF(${PACKAGE_NAME}_ENABLE_Ifpack)

#  PACKAGE_ADD_EXECUTABLE_AND_TEST(
#    tfqmr_hb
#    SOURCES test_tfqmr_hb.cpp createEpetraProblem.cpp
#    COMM serial mpi
#    ARGS
#      "--verbose --filename=orsirr1.hb"
#    STANDARD_PASS_OUTPUT 
#    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyTestTFQMRFiles
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/example/BlockGmres
    SOURCE_FILES orsirr1.hb orsirr1_scaled.hb
    EXEDEPS tfqmr_hb
    )

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
