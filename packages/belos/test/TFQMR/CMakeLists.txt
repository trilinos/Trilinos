
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
    
  ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)
  IF(Teuchos_ENABLE_COMPLEX)

    INCLUDE_DIRECTORIES(../MVOPTester)

    PACKAGE_ADD_EXECUTABLE(
      tfqmr_complex_hb
      SOURCES test_tfqmr_complex_hb.cpp 
      )

    PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyTestTFQMRComplexFiles
      SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/test/MVOPTester
      SOURCE_FILES mhd1280b.cua
      EXEDEPS tfqmr_complex_hb
      )

  ENDIF(Teuchos_ENABLE_COMPLEX)

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
