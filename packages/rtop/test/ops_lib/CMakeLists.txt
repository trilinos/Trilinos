INCLUDE(PackageAddExecutableAndTest)


INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/test/support)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
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
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  ARGS "--error-tol-slack=1e+2"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )
