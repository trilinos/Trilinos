IF (TrilinosCouplings_ENABLE_Tpetra AND Tpetra_Have_Kokkos_Refactor AND Trilinos_ENABLE_CXX11)

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

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
      Stokhos_ENABLE_TpetraKernels AND
      Stokhos_ENABLE_Tpetra )

    # The ensemble version requires the PCE scalar type too since it builds
    # a PCE through ensembles.
    IF(Stokhos_ENABLE_Ensemble_Scalar_Type AND
       Stokhos_ENABLE_PCE_Scalar_Type)
      TRIBITS_ADD_EXECUTABLE(
        fenl_ensemble
        SOURCES  fenl_ensemble.hpp SampleGrouping.hpp main_ensemble.cpp
        COMM serial mpi
        TESTONLYLIBS fenl
        )
    ENDIF()

    IF(Stokhos_ENABLE_PCE_Scalar_Type)
      TRIBITS_ADD_EXECUTABLE(
        fenl_pce
        SOURCES fenl_functors_pce.hpp fenl_pce.hpp main_pce.cpp
        COMM serial mpi
        TESTONLYLIBS fenl
        )
    ENDIF()

  ENDIF()

ENDIF ()
