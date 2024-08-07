

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../MVOPTester)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
    
  ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)
  IF(Teuchos_ENABLE_COMPLEX)

    TRIBITS_INCLUDE_DIRECTORIES(../MVOPTester)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      bicgstab_complex_hb
      SOURCES test_bicgstab_complex_hb.cpp 
      ARGS
        "--verbose --filename=mhd1280b.cua"
      )

    ASSERT_DEFINED(Anasazi_SOURCE_DIR)
    TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestBiCGStabComplexFiles
      SOURCE_DIR ${Anasazi_SOURCE_DIR}/testmatrices
      SOURCE_FILES mhd1280b.cua
      EXEDEPS bicgstab_complex_hb
    )

  ENDIF(Teuchos_ENABLE_COMPLEX)

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
