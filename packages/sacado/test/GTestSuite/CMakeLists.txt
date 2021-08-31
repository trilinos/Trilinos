INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/Fad)
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/utils)

SET(${PACKAGE_NAME}_TESTSUITE_SOURCES "")

SET(${PACKAGE_NAME}_TESTSUITE_SOURCES ${${PACKAGE_NAME}_TESTSUITE_SOURCES}
  FadUnitTests.hpp
  FadUnitTests.cpp
  CacheFadUnitTests.cpp
  ELRFadUnitTests.cpp
  ELRCacheFadUnitTests.cpp
  ComplexFadUnitTests.cpp
  ComplexCacheFadUnitTests.cpp
  ComplexELRFadUnitTests.cpp
  ComplexELRCacheFadUnitTests.cpp
  RealFadUnitTests.cpp
  RealCacheFadUnitTests.cpp
  RealELRFadUnitTests.cpp
  RealELRCacheFadUnitTests.cpp
  NestedFadUnitTests.hpp
  NestedFadUnitTests.cpp
  TraitsTests.hpp
  TraitsTests.cpp
  LogicalSparseUnitTests.cpp
  TayUnitTests.hpp
  TayUnitTests.cpp
  GTestUtils.hpp
  GTestSuite.cpp
  )

IF (Sacado_ENABLE_TeuchosNumerics)
  SET(${PACKAGE_NAME}_TESTSUITE_SOURCES ${${PACKAGE_NAME}_TESTSUITE_SOURCES}
     FadBLASUnitTests.hpp
     RealDFadBLASUnitTests.cpp
     RealDVFadBLASUnitTests.cpp
     ComplexDFadBLASUnitTests.cpp
     ComplexDVFadBLASUnitTests.cpp
     RealDFadLAPACKUnitTests.cpp
     )
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  GTestSuite
  SOURCES ${${PACKAGE_NAME}_TESTSUITE_SOURCES}
  ARGS
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "PASSED"
  )
