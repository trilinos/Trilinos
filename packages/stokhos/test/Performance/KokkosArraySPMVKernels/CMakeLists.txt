ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

IF (Stokhos_ENABLE_EpetraExt AND Stokhos_ENABLE_Kokkos)

  SET(CUDA_SOURCES "")
  IF (Stokhos_ENABLE_CUDA AND Kokkos_ENABLE_CUDA AND Stokhos_ENABLE_CUSPARSE)
    SET(CUDA_SOURCES TestCuda.cpp)
  ELSE()
    SET(CUDA_SOURCES TestCudaStub.cpp)
  ENDIF()
  TRIBITS_ADD_EXECUTABLE(
    SPMVTest
    SOURCES TestStochastic.hpp
            TestHost.cpp
            main.cpp
            ${CUDA_SOURCES}
    COMM serial mpi
    )

  IF(Kokkos_ENABLE_THREADS)
    TRIBITS_ADD_EXECUTABLE(
      SPMVHostScaling
      SOURCES TestStochastic.hpp
              HostScaling.cpp
      COMM serial mpi
      )
  ENDIF()

  TRIBITS_ADD_EXECUTABLE(
    SPMVEpetra
    SOURCES TestEpetra.cpp
    COMM serial mpi
    )

ENDIF()
