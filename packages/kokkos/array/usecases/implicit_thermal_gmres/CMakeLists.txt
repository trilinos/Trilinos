INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(TribitsAddExecutable)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_IMPLICIT_THERMAL_GMRES_HOST
  IMPLICIT_THERMAL_GMRES_HOST
  "Enable the host test for the Implicit Thermal GMRES usecase."
  OFF
  )

IF ( Kokkos_ENABLE_IMPLICIT_THERMAL_GMRES_HOST)
  MESSAGE(STATUS "Including Host test in Implicit Thermal GMRES usecase.")
  TRIBITS_ADD_EXECUTABLE(
    host_mini_test
    SOURCES mainhost.cpp testHost.cpp
    DEPLIBS kokkosdevicehost kokkosimpl 
  )
ENDIF (Kokkos_ENABLE_IMPLICIT_THERMAL_GMRES_HOST)

TRIBITS_ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_IMPLICIT_THERMAL_GMRES_PTHREAD
  IMPLICIT_THERMAL_GMRES_PTHREAD
  "Enable the pthread test for the Implicit Thermal GMRES usecase."
  OFF
  )

IF (Kokkos_ENABLE_IMPLICIT_THERMAL_GMRES_PTHREAD)
  MESSAGE(STATUS "Including Pthread test in Implicit Thermal GMRES usecase.")
  TRIBITS_ADD_EXECUTABLE(
    pthread_mini_test
    SOURCES mainpthread.cpp testPthread.cpp
    DEPLIBS kokkosdevicepthread kokkosdevicehost kokkosimpl
  )
ENDIF (Kokkos_ENABLE_IMPLICIT_THERMAL_GMRES_PTHREAD)

TRIBITS_ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_IMPLICIT_THERMAL_GMRES_TPI
  IMPLICIT_THERMAL_GMRES_TPI
  "Enable the tpi test for the Implicit Thermal GMRES usecase."
  OFF
  )
IF ( Kokkos_ENABLE_IMPLICIT_THERMAL_GMRES_TPI)
  MESSAGE(STATUS "Including TPI test in Implicit Thermal GMRES usecase.")
  
  TRIBITS_ADD_EXECUTABLE(
    tpi_mini_test
    SOURCES maintpi.cpp testTPI.cpp
    DEPLIBS kokkosdevicetpi kokkosdevicehost kokkosimpl
  )
ENDIF (Kokkos_ENABLE_IMPLICIT_THERMAL_GMRES_TPI)


TRIBITS_ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_IMPLICIT_THERMAL_GMRES_CUDA
  IMPLICIT_THERMAL_GMRES_CUDA
  "Enable the cuda test for the Implicit Thermal GMRES usecase."
  OFF
  )

IF ( Kokkos_ENABLE_IMPLICIT_THERMAL_GMRES_CUDA)
  MESSAGE(STATUS "Including Cuda test in Implicit Thermal GMRES usecase.")
  CUDA_COMPILE(CUDA_FILES testCuda.cu)
  TRIBITS_ADD_EXECUTABLE(
    cuda_mini_test
    SOURCES maincuda.cpp ${CUDA_FILES}
    DEPLIBS kokkosdevicecuda kokkosimpl
  )
ENDIF (Kokkos_ENABLE_IMPLICIT_THERMAL_GMRES_CUDA)
