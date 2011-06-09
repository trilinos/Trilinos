
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

APPEND_SET(TPI_HEADERS
    Kokkos_DeviceTPI_MDArrayView.hpp
    Kokkos_DeviceTPI_ParallelFor.hpp
    Kokkos_DeviceTPI_ParallelReduce.hpp
    Kokkos_DeviceTPI_ValueView.hpp
  )

APPEND_SET(TPI_SOURCES
    Kokkos_DeviceTPI.cpp
  )

PACKAGE_ADD_LIBRARY(
  kokkosdevicetpi
  HEADERS ${TPI_HEADERS}
  SOURCES ${TPI_SOURCES}
  DEPLIBS kokkos 
  )


