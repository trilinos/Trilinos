IF (TrilinosCouplings_ENABLE_Tpetra)

  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

  TRIBITS_ADD_LIBRARY(
    fenl
    HEADERS HexElement.hpp
            BoxElemFixture.hpp
            BoxElemPart.hpp
            fenl.hpp
            fenl_impl.hpp
            fenl_utils.hpp
            fenl_functors.hpp
            CGSolve.hpp
            BelosSolve.hpp
    SOURCES BoxElemPart.cpp
            fenl_utils.cpp
    TESTONLY
    )

  TRIBITS_ADD_EXECUTABLE(
    fenl
    SOURCES main.cpp
    COMM serial mpi
    TESTONLYLIBS fenl
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(copyFenlFiles
    DEST_FILES fenl.xml
    )

  IF( TrilinosCouplings_ENABLE_Stokhos AND
      Stokhos_ENABLE_Sacado AND
      Stokhos_ENABLE_KokkosKernels AND
      Stokhos_ENABLE_Tpetra )

    # The ensemble version requires the PCE scalar type too since it builds
    # a PCE through ensembles.
    IF(Stokhos_ENABLE_Ensemble_Scalar_Type AND
       Stokhos_ENABLE_PCE_Scalar_Type)
      TRIBITS_ADD_EXECUTABLE(
        fenl_ensemble
        SOURCES  fenl_ensemble.hpp SampleGrouping.hpp main_ensemble.cpp
                 VPS.hpp VPS.cpp VPS_ensemble.hpp VPS_ensemble.cpp
        COMM serial mpi
        TESTONLYLIBS fenl
        )

      IF(TrilinosCouplings_ENABLE_Belos AND TrilinosCouplings_ENABLE_MueLu)
        TRIBITS_ADD_TEST(
          fenl_ensemble
          NAME "fenl_ensemble_8x8x8_belos_muelu"
          COMM serial mpi
          #NUM_MPI_PROCS 4
          ARGS "--fixture=8x8x8 --belos --muelu --ensemble=${Stokhos_ETI_Default_Ensemble_Size} --print-its --unit-test --test-mean=2.0599472 --test-std-dev=8.9771555e-3"
          PASS_REGULAR_EXPRESSION "Test Passed"
          )
      ENDIF()

    ENDIF()

    IF(Stokhos_ENABLE_PCE_Scalar_Type)
      TRIBITS_ADD_EXECUTABLE(
        fenl_pce
        SOURCES fenl_functors_pce.hpp fenl_pce.hpp main_pce.cpp
        COMM serial mpi
        TESTONLYLIBS fenl
        )

      IF(TrilinosCouplings_ENABLE_Belos AND TrilinosCouplings_ENABLE_MueLu)
        TRIBITS_ADD_TEST(
          fenl_pce
          NAME "fenl_pce_8x8x8_belos_muelu_mean_based"
          COMM serial mpi
          #NUM_MPI_PROCS 4
          ARGS "--fixture=8x8x8 --belos --muelu --mean-based --print-its --unit-test --test-mean=2.0599472 --test-std-dev=8.9771555e-3"
          PASS_REGULAR_EXPRESSION "Test Passed"
          )
      ENDIF()
    ENDIF()

  ENDIF()

ENDIF ()
