
INCLUDE(TribitsAddExecutableAndTest)

ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

IF (Stokhos_ENABLE_KokkosCore)

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
  ELSE()
    SET(CUDA_SOURCES TestCudaStub.cpp)
  ENDIF()
  TRIBITS_ADD_EXECUTABLE(
    SPMVTestNew
    SOURCES TestStochastic.hpp
            TestHost.cpp
	    main.cpp
            ${CUDA_SOURCES}
    COMM serial mpi
    DEPLIBS ${DEPLIBS}
    )

  TRIBITS_ADD_EXECUTABLE(
    SPMVHostScalingNew
    SOURCES TestStochastic.hpp
            HostScaling.cpp
    COMM serial mpi
    DEPLIBS stokhos kokkoscore
    )

ENDIF()
