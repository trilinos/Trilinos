INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/Fad)

SET(${PACKAGE_NAME}_TESTSUITE_SOURCES "")

IF (TPL_ENABLE_ADOLC)

  SET(${PACKAGE_NAME}_TESTSUITE_SOURCES ${${PACKAGE_NAME}_TESTSUITE_SOURCES}
    TayUnitTests.hpp
    TayUnitTests.cpp
    )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Stokhos)

  SET(${PACKAGE_NAME}_TESTSUITE_SOURCES ${${PACKAGE_NAME}_TESTSUITE_SOURCES}
    HermiteUnitTests.hpp
    HermiteUnitTests.cpp
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
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "OK"
  )
