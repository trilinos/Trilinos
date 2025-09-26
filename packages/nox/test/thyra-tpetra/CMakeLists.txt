
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

ASSERT_DEFINED(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_THYRA
  NOX_ENABLE_Tpetra
  NOX_ENABLE_Stratimikos
  Stratimikos_ENABLE_Belos)

SET(UNIT_TEST_DRIVER UnitTestMainWithKokkos.cpp)

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_THYRA AND
    NOX_ENABLE_THYRA_TPETRA_ADAPTERS AND
    NOX_ENABLE_Stratimikos AND
    Stratimikos_ENABLE_Belos)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Thyra_Tpetra_VectorOpsTests
    SOURCES  ${UNIT_TEST_DRIVER} Tpetra_VectorOpsTests.cpp
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Thyra_Tpetra_MultiVectorOpsTests
    SOURCES  ${UNIT_TEST_DRIVER} Tpetra_MultiVectorOpsTests.cpp
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Thyra_Tpetra_OperatorTests
    SOURCES  ${UNIT_TEST_DRIVER} Tpetra_OperatorTests.cpp
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Thyra_Tpetra_1DFEM_ME_UnitTests
    SOURCES  ${UNIT_TEST_DRIVER} ME_Tpetra_1DFEM.hpp ME_Tpetra_1DFEM_def.hpp tTpetra_1DFEM_ME_UnitTests.cpp
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Thyra_Tpetra_1DFEM
    SOURCES  ${UNIT_TEST_DRIVER} ME_Tpetra_1DFEM.hpp ME_Tpetra_1DFEM_def.hpp Tpetra_1DFEM.cpp
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Thyra_Tpetra_Heq
    SOURCES  ${UNIT_TEST_DRIVER} ME_Tpetra_Heq.hpp ME_Tpetra_Heq_def.hpp Tpetra_Heq.cpp
    )

  IF(NOX_ENABLE_LOCA)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Thyra_Tpetra_HouseholderBorderedSolve
      SOURCES  ${UNIT_TEST_DRIVER} ME_Tpetra_1DFEM.hpp ME_Tpetra_1DFEM_def.hpp tTpetra_HouseholderBorderedSolve.cpp
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Thyra_Tpetra_HouseholderBorderedSolve_Continuation
      SOURCES  ${UNIT_TEST_DRIVER} ME_Tpetra_1DFEM.hpp ME_Tpetra_1DFEM_def.hpp tTpetra_HouseholderBorderedSolve_Continuation.cpp
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Thyra_Tpetra_HouseholderBorderedSolve_WithUVInPrec
      SOURCES  ${UNIT_TEST_DRIVER} ME_Tpetra_1DFEM.hpp ME_Tpetra_1DFEM_def.hpp tTpetra_HouseholderBorderedSolve_WithUVInPrec.cpp
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Thyra_Tpetra_ConstraintModelEvaluator
      SOURCES  ${UNIT_TEST_DRIVER} ME_Tpetra_1DFEM.hpp ME_Tpetra_1DFEM_def.hpp tTpetra_ConstraintModelEvaluator.cpp
      )

  ENDIF()

ENDIF()
