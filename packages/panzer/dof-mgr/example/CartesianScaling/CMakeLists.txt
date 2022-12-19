
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/cartesian_topology)

SET(ScalingTest_SOURCES
  main.cpp
  ${PACKAGE_SOURCE_DIR}/test/cartesian_topology/CartesianConnManager.cpp
  )

TRIBITS_ADD_EXECUTABLE(
  scaling
  SOURCES ${ScalingTest_SOURCES}
  )

TRIBITS_ADD_ADVANCED_TEST(
  scaling_test
  TEST_0 EXEC scaling
    PASS_REGULAR_EXPRESSION "Scaling test completed."
    XHOST trilinos-test.sandia.gov trilinos-test2.sandia.gov zan.sandia.gov
  )
