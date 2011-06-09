PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINRARY_DIR})

SET(HEADERS ${HEADERS}
 ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
    Kokkos_ArrayForwardDeclerations.hpp
    Kokkos_DeviceClear_macros.hpp
    Kokkos_DeviceCuda_macros.hpp
    Kokkos_DeviceCuda.hpp
    Kokkos_DeviceHost_macros.hpp
    Kokkos_DeviceHost.hpp
    Kokkos_DeviceTPI_macros.hpp
    Kokkos_DeviceTPI.hpp
    Kokkos_MDArrayView.hpp
    Kokkos_MemoryView.hpp
    Kokkos_MultiVectorView.hpp
    Kokkos_ParallelFor.hpp
    Kokkos_ParallelReduce.hpp
    Kokkos_ValueView.hpp
  )


INCLUDE(PackageLibraryMacros)
INCLUDE(AddOptionAndDefine)

ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_DeviceCuda
  HAVE_KOKKOS_DEVICECUDA
  "Enable arrays using CUDA Device."
  OFF
  )

ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_DeviceHost
  HAVE_KOKKOS_DEVICEHOST
  "Enable arrays using Host Device"
  OFF
  )

ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_DeviceTPI
  HAVE_KOKKOS_DEVICETPI
  "Enable arrays using TPI Device."
  OFF
  )

IF(HAVE_KOKKOS_DEVICECUDA)
  ADD_SUBDIRECTORY(DeviceCuda)
ENDIF()

IF(HAVE_KOKKOS_DEVICEHOST)
  ADD_SUBDIRECTORY(DeviceHost)
ENDIF()

IF(HAVE_KOKKOS_DEVICETPI)
  ADD_SUBDIRECTORY(DeviceTPI)
ENDIF()

ADD_SUBDIRECTORY(impl)


# No longer needed; Trilinos does this at the top level.
# IF (Kokkos_ENABLE_MDArray_Cuda)
#   FIND_PACKAGE(CUDA)
# ENDIF()
