INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/Fad)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../utils)

SET(${PACKAGE_NAME}_TESTSUITE_SOURCES "")

IF (TPL_ENABLE_ADOLC)

  SET(${PACKAGE_NAME}_TESTSUITE_SOURCES ${${PACKAGE_NAME}_TESTSUITE_SOURCES}
    CacheTaylorUnitTests.hpp
    CacheTaylorUnitTests.cpp
    TaylorUnitTests.hpp
    TaylorUnitTests.cpp
    )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Stokhos)

  SET(${PACKAGE_NAME}_TESTSUITE_SOURCES ${${PACKAGE_NAME}_TESTSUITE_SOURCES}
    HermiteUnitTests.hpp
    HermiteUnitTests.cpp
    )

ENDIF()


SET(${PACKAGE_NAME}_TESTSUITE_SOURCES ${${PACKAGE_NAME}_TESTSUITE_SOURCES}
  DFadUnitTests.hpp
  DFadUnitTests.cpp
  CacheDFadUnitTests.hpp
  CacheDFadUnitTests.cpp
  ELRDFadUnitTests.hpp
  ELRDFadUnitTests.cpp
  SFadUnitTests.hpp
  SFadUnitTests.cpp
  ELRSFadUnitTests.hpp
  ELRSFadUnitTests.cpp
  LogicalSparseUnitTests.hpp
  LogicalSparseUnitTests.cpp
  TestSuite.cpp
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  TestSuite
  SOURCES ${${PACKAGE_NAME}_TESTSUITE_SOURCES}
  DEPLIBS sacadoutils
  ARGS 
  COMM serial mpi
  )
