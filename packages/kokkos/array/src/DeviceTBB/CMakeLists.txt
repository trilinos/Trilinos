
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

APPEND_SET(TBB_HEADERS
    Kokkos_DeviceTBB_MDArrayView.hpp
    Kokkos_DeviceTBB_ParallelFor.hpp
    Kokkos_DeviceTBB_ParallelReduce.hpp
    Kokkos_DeviceTBB_ValueView.hpp
  )

APPEND_SET(TBB_SOURCES
    Kokkos_DeviceTBB.cpp
  )

PACKAGE_ADD_LIBRARY(
  kokkosdevicetbb
  HEADERS ${TBB_HEADERS}
  SOURCES ${TBB_SOURCES}
  DEPLIBS kokkos 
  )


