INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  projected_least_squares_solver
  SOURCES test_belos_projected_least_squares_solver.cpp
  ARGS "--verbose"
  )
