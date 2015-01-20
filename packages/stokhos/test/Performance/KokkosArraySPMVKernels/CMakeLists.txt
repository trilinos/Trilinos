ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

IF (Stokhos_ENABLE_EpetraExt)

  SET(CUDA_SOURCES "")
  IF (Stokhos_ENABLE_CUDA AND Kokkos_ENABLE_Cuda)
    SET(CUDA_SOURCES TestCuda.cpp)
  ENDIF()
  TRIBITS_ADD_EXECUTABLE(
    SPMVTest
    SOURCES TestStochastic.hpp
            TestHost.cpp
            main.cpp
            ${CUDA_SOURCES}
    COMM serial mpi
    )

  IF(Kokkos_ENABLE_Pthread)
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
