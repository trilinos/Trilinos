TRIBITS_ADD_EXECUTABLE_AND_TEST(
  hermite_example
  SOURCES simple/hermite_example.cpp
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "Example Passed"
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE(
  exp_moment_example
  SOURCES simple/exp_moment_example.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  recurrence_example
  SOURCES simple/recurrence_example.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  random_field_example
  SOURCES simple/random_field_example.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  basis_example
  SOURCES simple/basis.cpp
  COMM serial mpi
  )


TRIBITS_ADD_EXECUTABLE(
  stieltjes_example
  SOURCES dimension_reduction/stieltjes_example.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  stieltjes_example2
  SOURCES dimension_reduction/stieltjes_example2.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  stieltjes_example3
  SOURCES dimension_reduction/stieltjes_example3.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  stieltjes_example4
  SOURCES dimension_reduction/stieltjes_example4.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  stieltjes_coupled_example
  SOURCES dimension_reduction/stieltjes_coupled.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  gram_schmidt_example
  SOURCES dimension_reduction/gram_schmidt_example.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  gram_schmidt_example2
  SOURCES dimension_reduction/gram_schmidt_example2.cpp
  COMM serial mpi
  )

IF(Stokhos_ENABLE_PCE_Scalar_Type)

  TRIBITS_ADD_EXECUTABLE(
    gram_schmidt_example3
    SOURCES dimension_reduction/gram_schmidt_example3.cpp
    COMM serial mpi
    )

  TRIBITS_ADD_EXECUTABLE(
    stieltjes_example5
    SOURCES dimension_reduction/stieltjes_example5.cpp
    COMM serial mpi
    )

ENDIF()


IF(Stokhos_ENABLE_PCE_Scalar_Type)

  TRIBITS_ADD_EXECUTABLE(
    uq_handbook_pce_example
    SOURCES uq_handbook/pce_example.cpp
    COMM serial mpi
    )


  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    sacado_example
    SOURCES sacado/sacado_example.cpp
    ARGS
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    division_example
    SOURCES sacado/division_example.cpp
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Example Passed"
    NUM_MPI_PROCS 1
    )

ENDIF()

IF(Stokhos_ENABLE_Ensemble_Scalar_Type)
  TRIBITS_ADD_EXECUTABLE(
    uq_handbook_ensemble_example
    SOURCES uq_handbook/ensemble_example.cpp
    COMM serial mpi
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    sacado_ensemble_example
    SOURCES sacado/sacado_ensemble_example.cpp
    ARGS
    COMM serial mpi
    NUM_MPI_PROCS 1
    )
ENDIF()

IF(Stokhos_ENABLE_Cusp)

  TRIBITS_ADD_EXECUTABLE(
    cusp_sa_blockcg_example
    SOURCES cusp/cusp_sa_blockcg.cpp
    COMM serial mpi
    )
ENDIF()

IF(Stokhos_ENABLE_Experimental)
  ADD_SUBDIRECTORY(experimental)
ENDIF()
