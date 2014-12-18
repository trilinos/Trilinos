
IF (TrilinosCouplings_ENABLE_Tpetra AND Tpetra_Have_Kokkos_Refactor)

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
            VectorImport.hpp
            CGSolve.hpp
            BelosSolve.hpp
    SOURCES BoxElemPart.cpp
            fenl_utils.cpp
    TESTONLY
    )

  SET(SOURCES fenl.cpp main.cpp)
  IF (Kokkos_ENABLE_Cuda)
    # add MPI include dir to nvcc
    IF( TPL_ENABLE_MPI )
      CUDA_INCLUDE_DIRECTORIES( ${MPI_BASE_DIR}/include )
    ENDIF()
    CUDA_INCLUDE_DIRECTORIES( ${Tpetra_INCLUDE_DIRS} )
    CUDA_INCLUDE_DIRECTORIES( ${KokkosCore_INCLUDE_DIRS} )
    CUDA_INCLUDE_DIRECTORIES( ${KokkosContainers_INCLUDE_DIRS} )
    CUDA_INCLUDE_DIRECTORIES( ${TpetraKernels_INCLUDE_DIRS} )
    CUDA_INCLUDE_DIRECTORIES( ${KokkosCompat_INCLUDE_DIRS} )
    IF (TrilinosCouplings_ENABLE_MueLu)
      CUDA_INCLUDE_DIRECTORIES( ${MueLu_INCLUDE_DIRS} )
    ENDIF()
    CUDA_COMPILE( OBJECTS_CUDA fenl.cpp.cu )
    LIST( APPEND SOURCES ${OBJECTS_CUDA} )
  ENDIF ()
  TRIBITS_ADD_EXECUTABLE(
    fenl
    SOURCES ${SOURCES}
    COMM serial mpi
    TESTONLYLIBS fenl
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(copyFenlFiles
    DEST_FILES muelu.xml belos.xml
    )

  IF( TrilinosCouplings_ENABLE_Stokhos AND
      Stokhos_ENABLE_Sacado AND
      Stokhos_ENABLE_KokkosCore AND
      Stokhos_ENABLE_TpetraKernels AND
      Stokhos_ENABLE_Tpetra )

    SET( SOURCES fenl_ensemble.hpp
                 fenl_ensemble_macros.hpp
                 fenl_ensemble_threads.cpp
                 fenl_ensemble_openmp.cpp)
    IF (Kokkos_ENABLE_Cuda)
      CUDA_COMPILE( OBJECTS_CUDA fenl_ensemble_cuda.cu main_ensemble.cpp.cu )
      LIST( APPEND SOURCES ${OBJECTS_CUDA} )
    ELSE()
      LIST( APPEND SOURCES main_ensemble.cpp )
    ENDIF()
    TRIBITS_ADD_EXECUTABLE(
      fenl_ensemble
      SOURCES ${SOURCES}
      COMM serial mpi
      TESTONLYLIBS fenl
      )

    SET( SOURCES fenl_functors_pce.hpp fenl_pce.cpp )
    IF (Kokkos_ENABLE_Cuda)
      CUDA_COMPILE( OBJECTS_CUDA fenl_pce.cpp.cu main_pce.cpp.cu )
      LIST( APPEND SOURCES ${OBJECTS_CUDA} )
    ELSE()
      LIST( APPEND SOURCES main_pce.cpp )
    ENDIF()
    TRIBITS_ADD_EXECUTABLE(
      fenl_pce
      SOURCES ${SOURCES}
      COMM serial mpi
      TESTONLYLIBS fenl
      )

  ENDIF()

ENDIF ()
