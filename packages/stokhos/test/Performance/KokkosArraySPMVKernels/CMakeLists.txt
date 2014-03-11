INCLUDE(TribitsAddExecutableAndTest)

ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

IF (Stokhos_ENABLE_KokkosCore AND Stokhos_ENABLE_EpetraExt)

  SET(DEPLIBS stokhos kokkoscore)
  SET(CUDA_SOURCES "")
  IF (Stokhos_ENABLE_CUDA AND Kokkos_ENABLE_Cuda)
    GET_PROPERTY(Stokhos_CUDA_INC_DIRS DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
    CUDA_INCLUDE_DIRECTORIES(${Stokhos_CUDA_INC_DIRS})
    # add MPI include dir to nvcc
    IF( TPL_ENABLE_MPI )
      CUDA_INCLUDE_DIRECTORIES( ${MPI_BASE_DIR}/include )
    ENDIF()
    CUDA_COMPILE(CUDA_SOURCES TestCuda.cu)
  ENDIF()
  TRIBITS_ADD_EXECUTABLE(
    SPMVTest
    SOURCES TestStochastic.hpp
            TestHost.cpp
            main.cpp
            ${CUDA_SOURCES}
    COMM serial mpi
    DEPLIBS ${DEPLIBS}
    )

  IF(Kokkos_ENABLE_Pthread)
    TRIBITS_ADD_EXECUTABLE(
      SPMVHostScaling
      SOURCES TestStochastic.hpp
              HostScaling.cpp
      COMM serial mpi
      DEPLIBS stokhos kokkoscore
      )
  ENDIF()

  TRIBITS_ADD_EXECUTABLE(
    SPMVEpetra
    SOURCES TestEpetra.cpp
    COMM serial mpi
    DEPLIBS stokhos
    )

ENDIF()
