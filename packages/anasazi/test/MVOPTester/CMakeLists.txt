INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
    
  ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)
  IF(Teuchos_ENABLE_COMPLEX)

    PACKAGE_ADD_EXECUTABLE(
      MVOPComplexTester
      SOURCES cxx_main_complex.cpp MyMultiVec.hpp MyOperator.hpp MyBetterOperator.hpp
      COMM serial mpi
      )

    PACKAGE_ADD_TEST(
      MVOPComplexTester
      COMM serial
      ARGS
        "--verbose --filename=mhd1280b.cua"
      STANDARD_PASS_OUTPUT
      )

    PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyTestMVOPTesterComplexFiles
      SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/anasazi/testmatrices
      SOURCE_FILES mhd1280b.cua 
      EXEDEPS MVOPComplexTester
      )

  ENDIF(Teuchos_ENABLE_COMPLEX)

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
