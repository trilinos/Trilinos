
INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddTest)
INCLUDE(PackageCopyFilesToBinaryDir)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
    
  ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)
  IF(Teuchos_ENABLE_COMPLEX)

    INCLUDE_DIRECTORIES(../MVOPTester)

    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      OrthoManagerComplexTester
      SOURCES cxx_main_complex.cpp 
      COMM serial
      ARGS
        "--verbose --ortho=SVQB"
	"--verbose --filename=mhd1280b.cua --ortho=SVQB"
        "--verbose --ortho=Basic"
	"--verbose --filename=mhd1280b.cua --ortho=Basic"
        "--verbose --ortho=ICGS"
	"--verbose --filename=mhd1280b.cua --ortho=ICGS"
      STANDARD_PASS_OUTPUT
      )

    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      OrthoManagerComplexMatTester
      SOURCES cxx_mattest_complex.cpp 
      COMM serial
      ARGS
        "--verbose --ortho=SVQB"
	"--verbose --filename=mhd1280b.cua --ortho=SVQB"
        "--verbose --ortho=Basic"
	"--verbose --filename=mhd1280b.cua --ortho=Basic"
        "--verbose --ortho=ICGS"
	"--verbose --filename=mhd1280b.cua --ortho=ICGS"
      STANDARD_PASS_OUTPUT
      )

    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      OrthoManagerComplexGenTester
      SOURCES cxx_gentest_complex.cpp 
      COMM serial
      ARGS
        "--verbose"
        "--verbose --filename=mhd1280b.cua"
      STANDARD_PASS_OUTPUT
      )

    PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyTestOrthoManagerTesterComplexFiles
      SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/anasazi/test/MVOPTester
      SOURCE_FILES mhd1280b.cua 
      EXEDEPS OrthoManagerComplexTester
      )

  ENDIF(Teuchos_ENABLE_COMPLEX)

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
