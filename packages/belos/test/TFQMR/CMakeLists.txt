

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../MVOPTester)

ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)
IF(Teuchos_ENABLE_COMPLEX)

  TRIBITS_INCLUDE_DIRECTORIES(../MVOPTester)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    tfqmr_complex_hb
    SOURCES test_tfqmr_complex_hb.cpp
    ARGS
      "--verbose --filename=mhd1280b.cua"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    pseudo_tfqmr_complex_hb
    SOURCES test_pseudo_tfqmr_complex_hb.cpp
    ARGS
      "--verbose --filename=mhd1280b.cua"
      "--verbose --num-rhs=5 --filename=mhd1280b.cua"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    tfqmr_complex_diag
    SOURCES test_tfqmr_complex_diag.cpp
    ARGS
      "--verbose"
      "--verbose --pseudo"
    COMM serial
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestTFQMRComplexFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
    SOURCE_FILES mhd1280b.cua
    EXEDEPS tfqmr_complex_hb
  )

ENDIF(Teuchos_ENABLE_COMPLEX)
