INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddTest)
INCLUDE(PackageCopyFilesToBinaryDir)

ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)

INCLUDE_DIRECTORIES(../MVOPTester)

IF(Teuchos_ENABLE_COMPLEX)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
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

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    IRTR_complex_test
    SOURCES cxx_main_complex.cpp
    COMM serial mpi
    ARGS
      "--sort=LR --verbose --hefty --filename=mhd1280b.cua"
      "--sort=LR --verbose --skinny --filename=mhd1280b.cua"
    STANDARD_PASS_OUTPUT
    )

    PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyTestIRTRComplexFiles
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/anasazi/testmatrices
    SOURCE_FILES mhd1280b.cua
    EXEDEPS IRTR_complex_test
    )

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)

ENDIF()
