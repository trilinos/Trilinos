
INCLUDE(PackageAddExecutableAndTest)

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
  ImplicitRK
  )

FOREACH(TEST_NAME ${TEST_NAMES})
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    ${TEST_NAME}_ConvergenceTest
    SOURCES
      ${COMMON_SOURCES}
      Rythmos_${TEST_NAME}_ConvergenceTest.cpp
    DEPLIBS rythmos_test_models
    COMM serial mpi  
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

ENDFOREACH()
