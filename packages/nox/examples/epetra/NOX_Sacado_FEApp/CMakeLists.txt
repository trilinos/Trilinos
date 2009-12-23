
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA AND NOX_ENABLE_Sacado AND
   NOX_ENABLE_EpetraExt AND Sacado_ENABLE_FEApp)
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Sacado_FEApp_Brusselator_Example
    SOURCES brusselator.C
    COMM serial mpi
    NUM_MPI_PROCS 2
    ARGS -v
    PASS_REGULAR_EXPRESSION "Test Passed"
    )

  IF(Sacado_ENABLE_Stokhos)
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      Sacado_FEApp_PCE_Example
      SOURCES pce.C ENAT_NOXSolver.hpp ENAT_NOXSolver.cpp
      COMM serial mpi
      NUM_MPI_PROCS 2
      ARGS
      PASS_REGULAR_EXPRESSION "Test Passed"
      )
    PACKAGE_ADD_EXECUTABLE(
      Sacado_FEApp_PCE_KL_Example
      SOURCES pce_KL.C ENAT_NOXSolver.hpp ENAT_NOXSolver.cpp
      COMM serial mpi
      )
    PACKAGE_ADD_EXECUTABLE(
      Jacobian_KL_Example
      SOURCES jacobian_KL.cpp ENAT_NOXSolver.hpp ENAT_NOXSolver.cpp
      COMM serial mpi
      )
  ENDIF()

ENDIF()
