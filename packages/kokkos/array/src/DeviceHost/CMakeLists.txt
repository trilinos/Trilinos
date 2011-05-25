
INCLUDE(PackageLibraryMacros)
INCLUDE(AddOptionAndDefine)

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINRARY_DIR})

SET(HEADERS ${HEADERS}
 ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
    Kokkos_DeviceHost_MDArrayView.hpp
    Kokkos_DeviceHost_MultiVectorView.hpp
    Kokkos_DeviceHost_ParallelFor.hpp
    Kokkos_DeviceHost_ParallelReduce.hpp
  )

APPEND_SET(SOURCES
    Kokkos_DeviceHost.cpp
  )

PACKAGE_ADD_LIBRARY(
  kokkosdevicehost
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS kokkos
  )


# No longer needed; Trilinos does this at the top level.
# IF (Kokkos_ENABLE_MDArray_Cuda)
#   FIND_PACKAGE(CUDA)
# ENDIF()
