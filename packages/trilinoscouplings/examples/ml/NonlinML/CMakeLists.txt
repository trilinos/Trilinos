
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ml_nox_1Delasticity_example
  SOURCES Basis.cpp
          ml_nox_1Delasticity_example.cpp
          FiniteElementProblem.cpp
          Problem_Interface.cpp
  ARGS
  COMM serial mpi
  )
