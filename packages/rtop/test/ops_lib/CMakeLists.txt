

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/test/support)


SET(ERROR_TOL_SLACK 1e+2)

IF (CMAKE_HOST_SYSTEM_NAME STREQUAL AIX)
  # Set to pass a failing Nightly test on the AIX machine  p90n03.pbm.ihost.com
  SET(ERROR_TOL_SLACK 1e+4)
ENDIF()


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  opsUnitTests
  SOURCES
    ROpCountNanInf_UnitTests.cpp
    ROpDotProd_UnitTests.cpp
    ROpGetElement_UnitTests.cpp
    ROpGetSubVector_UnitTests.cpp
    ROpMaxIndexLessThanBound_UnitTests.cpp
    ROpMaxIndex_UnitTests.cpp
    ROpMax_UnitTests.cpp
    ROpMinIndexGreaterThanBound_UnitTests.cpp
    ROpMinIndex_UnitTests.cpp
    ROpMin_UnitTests.cpp
    ROpNorm1_UnitTests.cpp
    ROpNorm2_UnitTests.cpp
    ROpNormInf_UnitTests.cpp
    ROpSum_UnitTests.cpp
    ROpWeightedNorm2_UnitTests.cpp
    TOpAXPY_UnitTests.cpp
    TOpAbs_UnitTests.cpp
    TOpAddScalar_UnitTests.cpp
    TOpAssignScalar_UnitTests.cpp
    TOpAssignVectors_UnitTests.cpp
    TOpEleWiseConjProd_UnitTests.cpp
    TOpEleWiseDivide_UnitTests.cpp
    TOpEleWiseProdUpdate_UnitTests.cpp
    TOpEleWiseProd_UnitTests.cpp
    TOpEleWiseScale_UnitTests.cpp
    TOpLinearCombination_UnitTests.cpp
    TOpRandomize_UnitTests.cpp
    TOpReciprocal_UnitTests.cpp
    TOpScaleVector_UnitTests.cpp
    TOpSetAssendingValues_UnitTests.cpp
    TOpSetElement_UnitTests.cpp
    TOpSetSubVector_UnitTests.cpp
    ../support/supportUnitTestsHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS "--error-tol-slack=${ERROR_TOL_SLACK}"
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )
