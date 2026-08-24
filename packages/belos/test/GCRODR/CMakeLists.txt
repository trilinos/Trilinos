

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../MVOPTester)

ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)
IF(Teuchos_ENABLE_COMPLEX)

  TRIBITS_INCLUDE_DIRECTORIES(../MVOPTester)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    gcrodr_complex_hb
    SOURCES test_gcrodr_complex_hb.cpp
    ARGS
      "--verbose --filename=mhd1280b.cua"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
      solution_projection_complex_hb
      SOURCES test_solution_projection_complex_hb.cpp
      ARGS
        "--verbose --filename=mhd1280b.cua --epsilon=1e-4 --projection-reduction-tol=1e-3 --require-projection-beats-recycled --min-iteration-improvement=200"
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestGCRODRComplexFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
    SOURCE_FILES mhd1280b.cua
    EXEDEPS gcrodr_complex_hb
  )

ENDIF(Teuchos_ENABLE_COMPLEX)
