
ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)
IF(Teuchos_ENABLE_COMPLEX)

  TRIBITS_ADD_EXECUTABLE(
    MVOPComplexTester
    SOURCES cxx_main_complex.cpp MyMultiVec.hpp MyOperator.hpp MyBetterOperator.hpp
    COMM serial mpi
    )

  TRIBITS_ADD_TEST(
    MVOPComplexTester
    COMM serial
    ARGS
      "--verbose --filename=mhd1280b.cua"
    STANDARD_PASS_OUTPUT
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestMVOPTesterComplexFiles
    SOURCE_DIR ${PACKAGE_SOURCE_DIR}/testmatrices
    SOURCE_FILES mhd1280b.cua
    EXEDEPS MVOPComplexTester
    )

ENDIF(Teuchos_ENABLE_COMPLEX)
