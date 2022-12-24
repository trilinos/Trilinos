
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
    
  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_COMPLEX)
  # FIXME (mfh 04 Aug 2010)
  IF ({${PACKAGE_NAME}_ENABLE_COMPLEX)

    INCLUDE_DIRECTORIES(../MVOPTester)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
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

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
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

    TRIBITS_ADD_EXECUTABLE(
      OrthoManagerComplexGenTester
      SOURCES cxx_gentest_complex.cpp 
      COMM serial mpi
      )

    TRIBITS_ADD_TEST(
      OrthoManagerComplexGenTester
      COMM serial
      ARGS
        "--verbose"
        "--verbose --filename=mhd1280b.cua"
      STANDARD_PASS_OUTPUT
      )

    TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestOrthoManagerTesterComplexFiles
      SOURCE_DIR ${PACKAGE_SOURCE_DIR}/testmatrices
      SOURCE_FILES mhd1280b.cua 
      EXEDEPS OrthoManagerComplexGenTester
      )

  ENDIF ()

ENDIF ()
