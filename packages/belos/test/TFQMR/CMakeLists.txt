
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../MVOPTester)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
    
  ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)
  IF(Teuchos_ENABLE_COMPLEX)

    TRIBITS_INCLUDE_DIRECTORIES(../MVOPTester)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      tfqmr_complex_hb
      SOURCES test_tfqmr_complex_hb.cpp 
      ARGS
        "--verbose --filename=mhd1280b.cua"
      )

    TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestTFQMRComplexFiles
      SOURCE_DIR ${Belos_SOURCE_DIR}/test/MVOPTester
      SOURCE_FILES mhd1280b.cua mhd1280a.cua
      EXEDEPS tfqmr_complex_hb
      )

  ENDIF(Teuchos_ENABLE_COMPLEX)

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
