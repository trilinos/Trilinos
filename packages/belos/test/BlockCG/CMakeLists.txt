
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../MVOPTester)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
    
  ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)
  IF(Teuchos_ENABLE_COMPLEX)

    INCLUDE_DIRECTORIES(../MVOPTester)

    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      bl_cg_complex_hb
      SOURCES test_bl_cg_complex_hb.cpp 
      ARGS
        "--verbose --filename=mhd1280b.cua"
      )

    PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyTestBlockCGComplexFiles
      SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/test/MVOPTester
      SOURCE_FILES mhd1280b.cua
      EXEDEPS bl_cg_complex_hb
      )

  ENDIF(Teuchos_ENABLE_COMPLEX)

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
