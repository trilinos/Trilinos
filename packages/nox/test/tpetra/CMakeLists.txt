
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

ASSERT_DEFINED(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_THYRA
  NOX_ENABLE_Tpetra
  NOX_ENABLE_Stratimikos
  Stratimikos_ENABLE_Belos)

SET(UNIT_TEST_DRIVER UnitTestMainWithKokkos.cpp)

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_THYRA AND
    NOX_ENABLE_THYRA_TPETRA_ADAPTERS AND
    NOX_ENABLE_Stratimikos AND
    Stratimikos_ENABLE_Belos)

  #TRIBITS_ADD_EXECUTABLE_AND_TEST(
  #  Tpetra_GroupTests
  #  SOURCES  ${UNIT_TEST_DRIVER} ME_Tpetra_1DFEM.hpp ME_Tpetra_1DFEM_def.hpp Tpetra_GroupTests.cpp
  #  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tpetra_VectorOpsTests
    SOURCES  ${UNIT_TEST_DRIVER} Tpetra_VectorOpsTests.cpp
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tpetra_MultiVectorOpsTests
    SOURCES  ${UNIT_TEST_DRIVER} Tpetra_MultiVectorOpsTests.cpp
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tpetra_OperatorTests
    SOURCES  ${UNIT_TEST_DRIVER} Tpetra_OperatorTests.cpp
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Belos_MultiVecTraitsTests
    SOURCES  ${UNIT_TEST_DRIVER} Belos_MultiVecTraitsTests.cpp
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Belos_OperatorTraitsTests
    SOURCES  ${UNIT_TEST_DRIVER} Belos_OperatorTraitsTests.cpp
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Belos_LOWSTests
    SOURCES  ${UNIT_TEST_DRIVER} Belos_LOWSTests.cpp
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tpetra_1DFEM
    SOURCES  ${UNIT_TEST_DRIVER} ME_Tpetra_1DFEM.hpp ME_Tpetra_1DFEM_def.hpp Tpetra_1DFEM.cpp
    )

ENDIF()
