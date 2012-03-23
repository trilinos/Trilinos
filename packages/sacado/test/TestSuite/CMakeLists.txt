INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/Fad)
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/utils)

SET(${PACKAGE_NAME}_TESTSUITE_SOURCES "")

IF (TPL_ENABLE_ADOLC)

  SET(${PACKAGE_NAME}_TESTSUITE_SOURCES ${${PACKAGE_NAME}_TESTSUITE_SOURCES}
    TayUnitTests.hpp
    TayUnitTests.cpp
    )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Teuchos)

  SET(${PACKAGE_NAME}_TESTSUITE_SOURCES ${${PACKAGE_NAME}_TESTSUITE_SOURCES}
    FadBLASUnitTests.hpp
    RealDFadBLASUnitTests.cpp
    RealDVFadBLASUnitTests.cpp
    ComplexDFadBLASUnitTests.cpp
    ComplexDVFadBLASUnitTests.cpp
    FadLAPACKUnitTests.hpp
    RealDFadLAPACKUnitTests.cpp
    )

ENDIF()


SET(${PACKAGE_NAME}_TESTSUITE_SOURCES ${${PACKAGE_NAME}_TESTSUITE_SOURCES}
  FadUnitTests.hpp
  FadUnitTests2.hpp
  FadUnitTests.cpp
  RealFadUnitTests.cpp
  ComplexFadUnitTests.cpp
  ELRFadUnitTests.cpp
  RealELRFadUnitTests.cpp
  ComplexELRFadUnitTests.cpp
  CacheFadUnitTests.cpp
  RealCacheFadUnitTests.cpp
  ComplexCacheFadUnitTests.cpp
  ELRCacheFadUnitTests.cpp
  RealELRCacheFadUnitTests.cpp
  ComplexELRCacheFadUnitTests.cpp
  LogicalSparseUnitTests.hpp
  LogicalSparseUnitTests.cpp
  TraitsTests.hpp
  TraitsTests.cpp
  TestSuite.cpp
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TestSuite
  SOURCES ${${PACKAGE_NAME}_TESTSUITE_SOURCES}
  ARGS 
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "OK"
  )
