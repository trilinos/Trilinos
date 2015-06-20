
ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ../MVOPTester)

IF(Teuchos_ENABLE_COMPLEX)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    IRTR_complex_lap_test
    SOURCES cxx_main_complex_lap.cpp
    ARGS 
      "--verbose --hefty --fakeM --sort=SR"
      "--verbose --skinny --noM --sort=SR"
      "--verbose --hefty --fakeM --sort=SR"
      "--verbose --skinny --noM --sort=SR"
      "--verbose --hefty --fakeM --sort=LR"
      "--verbose --skinny --noM --sort=LR" 
      "--verbose --hefty --fakeM --sort=LR"
      "--verbose --skinny --noM --sort=LR"
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    IRTR_complex_test
    SOURCES cxx_main_complex.cpp
    COMM serial mpi
    ARGS
      "--sort=LR --verbose --hefty --filename=mhd1280b.cua"
      "--sort=LR --verbose --skinny --filename=mhd1280b.cua"
    STANDARD_PASS_OUTPUT
    )

    TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestIRTRComplexFiles
    SOURCE_DIR ${PACKAGE_SOURCE_DIR}/testmatrices
    SOURCE_FILES mhd1280b.cua
    EXEDEPS IRTR_complex_test
    )

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)

ENDIF()
