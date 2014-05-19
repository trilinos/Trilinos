INCLUDE(TribitsAddExecutableAndTest)

IF (TrilinosCouplings_ENABLE_Tpetra AND Tpetra_Have_Kokkos_Refactor)

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

  SET(SOURCES BoxElemPart.cpp fenl.cpp main.cpp)

  IF (Kokkos_ENABLE_Cuda)

    # add MPI include dir to nvcc
    IF( TPL_ENABLE_MPI )
      CUDA_INCLUDE_DIRECTORIES( ${MPI_BASE_DIR}/include )
    ENDIF()

    CUDA_INCLUDE_DIRECTORIES( ${Tpetra_INCLUDE_DIRS} )
    CUDA_INCLUDE_DIRECTORIES( ${KokkosCore_INCLUDE_DIRS} )
    CUDA_INCLUDE_DIRECTORIES( ${KokkosContainers_INCLUDE_DIRS} )
    CUDA_INCLUDE_DIRECTORIES( ${KokkosLinAlg_INCLUDE_DIRS} )
    CUDA_INCLUDE_DIRECTORIES( ${KokkosCompat_INCLUDE_DIRS} )
    IF (TrilinosCouplings_ENABLE_MueLu)
      CUDA_INCLUDE_DIRECTORIES( ${MueLu_INCLUDE_DIRS} )
    ENDIF()

    CUDA_COMPILE( OBJECTS_CUDA fenl.cpp.cu )

    LIST( APPEND SOURCES ${OBJECTS_CUDA} )

  ENDIF () # Kokkos_ENABLE_Cuda

  TRIBITS_ADD_EXECUTABLE(
    fenl
    SOURCES ${SOURCES}
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(copyFenlFiles
    DEST_FILES muelu.xml
    )

  IF( TrilinosCouplings_ENABLE_Stokhos AND
      Stokhos_ENABLE_Sacado AND
      Stokhos_ENABLE_KokkosCore AND
      Stokhos_ENABLE_KokkosLinAlg )

    SET( SOURCES BoxElemPart.cpp fenl_ensemble.cpp )

    IF (Kokkos_ENABLE_Cuda)
      CUDA_COMPILE( OBJECTS_CUDA fenl_ensemble.cpp.cu main_ensemble.cpp.cu )
      LIST( APPEND SOURCES ${OBJECTS_CUDA} )
    ELSE()
      LIST( APPEND SOURCES main_ensemble.cpp )
    ENDIF()

    TRIBITS_ADD_EXECUTABLE(
      fenl_ensemble
      SOURCES ${SOURCES}
      COMM serial mpi
      DEPLIBS kokkoscorelinalg kokkoscontainers kokkoscore
      )

  ENDIF()

ENDIF ()
