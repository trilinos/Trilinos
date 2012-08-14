
INCLUDE(TribitsAddExecutableAndTest)

ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

SET(COMMON_SOURCES
  Rythmos_ConvergenceTestHelpers.cpp
  Rythmos_ConvergenceTest.cpp
  )

SET(TEST_NAMES
  BackwardEuler
  ForwardEuler
  ImplicitBDF
  ExplicitRK
  IntegratorBuilder
  )

IF (${PACKAGE_NAME}_ENABLE_Stratimikos  AND ${PACKAGE_NAME}_ENABLE_Belos)
#diabling this test for trilinos-test2 because of bug #5485
#I had to change the logic here a bit to disable just the single
#test that I wanted to on trilinos-test2. Instead of just adding
#the test name to a list and then looping over it instead we add this test
#separately. If anything changes in the TRIBITS_ADD_EXECUTABLE_AND_TEST
#inside the foreach loop then it needs to be made here as well. Hopefully
#this will be a short lived change.
#bmp

#  SET(TEST_NAMES ${TEST_NAMES}  ImplicitRK )
  SET(TEST_NAME ImplicitRK)
  
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ${TEST_NAME}_ConvergenceTest
    SOURCES
      ${COMMON_SOURCES}
      Rythmos_${TEST_NAME}_ConvergenceTest.cpp
    DEPLIBS rythmos_test_models
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    XHOST trilinos-test2.sandia.gov zan.sandia.gov
    )
ENDIF()

FOREACH(TEST_NAME ${TEST_NAMES})
  
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ${TEST_NAME}_ConvergenceTest
    SOURCES
      ${COMMON_SOURCES}
      Rythmos_${TEST_NAME}_ConvergenceTest.cpp
    DEPLIBS rythmos_test_models
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

ENDFOREACH()

