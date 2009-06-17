
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA AND NOX_ENABLE_Sacado AND
   NOX_ENABLE_EpetraExt AND Sacado_ENABLE_FEApp)
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    NOX_Sacado_FEApp_Brusselator_Example
    SOURCES brusselator.C
    COMM serial mpi
    NUM_MPI_PROCS 2
    ARGS -v
    PASS_REGULAR_EXPRESSION "Test Passed"
    )

  IF(Sacado_ENABLE_Stokhos AND Stokhos_ENABLE_Quadrature)
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      NOX_Sacado_FEApp_PCE_Example
      SOURCES pce.C FEApp_NOXSolver.hpp FEApp_NOXSolver.cpp
      COMM serial mpi
      NUM_MPI_PROCS 2
      ARGS
      PASS_REGULAR_EXPRESSION "Test Passed"
      )
  ENDIF()

ENDIF()
