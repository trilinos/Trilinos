
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../MVOPTester)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
    
  ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)
  IF(Teuchos_ENABLE_COMPLEX)

    TRIBITS_INCLUDE_DIRECTORIES(../MVOPTester)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      bl_cg_complex_hb
      SOURCES test_bl_cg_complex_hb.cpp 
      ARGS
        "--verbose --filename=mhd1280b.cua"
      XHOST KALLIKRATES
      )

    IF(NOT ${${PROJECT_NAME}_HOSTNAME} STREQUAL "KALLIKRATES")
      TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestBlockCGComplexFiles
        SOURCE_DIR ${Belos_SOURCE_DIR}/test/MVOPTester
        SOURCE_FILES mhd1280b.cua
        EXEDEPS bl_cg_complex_hb
        )
    ENDIF()

  ENDIF(Teuchos_ENABLE_COMPLEX)

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
