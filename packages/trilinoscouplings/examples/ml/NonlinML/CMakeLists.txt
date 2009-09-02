
INCLUDE(PackageAddExecutableAndTest)

IF(HAVE_TRILINOSCOUPLINGS_ML AND HAVE_TRILINOSCOUPLINGS_NOX)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    ml_nox_1Delasticity_example
    SOURCES Basis.cpp
            ml_nox_1Delasticity_example.cpp
            FiniteElementProblem.cpp
            Problem_Interface.cpp
    ARGS
    COMM serial mpi
    )
ENDIF()
