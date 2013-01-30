
INCLUDE(TribitsAddExecutableAndTest)

ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

IF (Stokhos_ENABLE_KokkosArray)

  SET(KOKKOS_ARRAY_UQ ${CMAKE_CURRENT_SOURCE_DIR}/../../../../kokkos/array/research/uq_on_gpgpu/src)
  INCLUDE_DIRECTORIES(${KOKKOS_ARRAY_UQ})
  SET(DEPLIBS stokhos kokkosarray_impl kokkosarray_devicehost)
  SET(CUDA_SOURCES "")
  IF (Stokhos_ENABLE_CUDA AND KokkosArray_ENABLE_DeviceCuda)
    GET_PROPERTY(Stokhos_CUDA_INC_DIRS DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
    CUDA_INCLUDE_DIRECTORIES(${Stokhos_CUDA_INC_DIRS})
    CUDA_COMPILE(CUDA_SOURCES TestCuda.cu)
    SET(DEPLIBS ${DEPLIBS} kokkosarray_devicecuda)
  ELSE()
    SET(CUDA_SOURCES TestCudaStub.cpp)
  ENDIF()
  TRIBITS_ADD_EXECUTABLE(
    SPMVTest
    SOURCES TestStochastic.hpp
            TestHost.cpp
	    main.cpp
            ${KOKKOS_ARRAY_UQ}/impl/KokkosArray_LegendrePolynomial.cpp
            ${KOKKOS_ARRAY_UQ}/impl/KokkosArray_ProductTensorLegendre.cpp
            ${CUDA_SOURCES}
    COMM serial mpi
    DEPLIBS ${DEPLIBS}
    )
ENDIF()