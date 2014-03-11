
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

#IF (${PACKAGE_NAME}_ENABLE_Stratimikos  AND ${PACKAGE_NAME}_ENABLE_Belos)
#  SET(TEST_NAMES ${TEST_NAMES}  ImplicitRK )
#ENDIF()

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

