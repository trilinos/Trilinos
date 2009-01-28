
INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddTest)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  MVOPTester
  SOURCES cxx_main.cpp 
  ARGS "-v"
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
    
  ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)
  IF(Teuchos_ENABLE_COMPLEX)

    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      MVOPComplexTester
      SOURCES cxx_main_complex.cpp MyMultiVec.hpp MyOperator.hpp MyBetterOperator.hpp
      COMM serial mpi
      ARGS
        "--verbose --filename=mhd1280b.cua"
      )

    PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyTestMVOPTesterComplexFiles
      SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/test/MVOPTester
      SOURCE_FILES mhd1280b.cua 
      EXEDEPS MVOPComplexTester
      )

  ENDIF(Teuchos_ENABLE_COMPLEX)

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
