INCLUDE(TribitsAddExecutableAndTest)

ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

IF (Stokhos_ENABLE_KokkosCore)

  SET(CUDA_SOURCES "")
  IF (Stokhos_ENABLE_CUDA AND Kokkos_ENABLE_Cuda)
    GET_PROPERTY(Stokhos_CUDA_INC_DIRS DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
    CUDA_INCLUDE_DIRECTORIES(${Stokhos_CUDA_INC_DIRS})
    # add MPI include dir to nvcc
    IF( TPL_ENABLE_MPI )
      CUDA_INCLUDE_DIRECTORIES( ${MPI_BASE_DIR}/include )
    ENDIF()
    CUDA_COMPILE(CUDA_SOURCES TestAssembly_Cuda.cu)
  ENDIF()
  TRIBITS_ADD_EXECUTABLE(
    MPVectorAssemblyTest
    SOURCES BoxElemFixture.hpp
            BoxElemPart.hpp
            BoxElemPart.cpp
            HexElement.hpp
            fenl_functors.hpp
            TestAssembly.hpp
            TestAssembly.cpp
            TestAssembly_Host.cpp
            ${CUDA_SOURCES}
    COMM serial mpi
    DEPLIBS stokhos kokkoscore
    )

ENDIF()
