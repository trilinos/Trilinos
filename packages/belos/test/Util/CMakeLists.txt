

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../MVOPTester)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
    
    TRIBITS_INCLUDE_DIRECTORIES(../MVOPTester)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Aztec_to_Belos_param
      SOURCES Aztec_to_Belos_param.cpp 
      ARGS
        " args ignored"
      )

#    ASSERT_DEFINED(Anasazi_SOURCE_DIR)
#    TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestBlockCGComplexFiles
#      SOURCE_DIR ${Anasazi_SOURCE_DIR}/testmatrices
#      SOURCE_FILES mhd1280b.cua 
#      EXEDEPS minres_complex_hb
#    )

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
