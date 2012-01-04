INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddOptionAndDefine)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS "")
SET(SOURCES "")

FILE(GLOB HEADERS *.hpp)

LIST(APPEND SOURCES Kokkos_DeviceNUMA.cpp)

#-----------------------------------------------------------------------------
# Hardware locality detection functionality depends upon TPL availability

IF( TPL_ENABLE_HWLOC )
#
#  http://www.open-mpi.org/projects/hwloc/
#  Version 1.3 or greater
#
LIST(APPEND SOURCES Kokkos_DeviceNUMA_hwloc.cpp)
#
ELSE()
LIST(APPEND SOURCES Kokkos_DeviceNUMA_hwloc_unavailable.cpp)
ENDIF()

#-----------------------------------------------------------------------------
# Thread spawning functionality depend upon TPL availability

IF( TPL_ENABLE_Pthread )
LIST(APPEND SOURCES Kokkos_DeviceNUMA_pthread.cpp)
ELSEIF( TPL_ENABLE_Winthread )
LIST(APPEND SOURCES Kokkos_DeviceNUMA_winthread.cpp)
ELSE()
LIST(APPEND SOURCES Kokkos_DeviceNUMA_nothread.cpp)
ENDIF()

#-----------------------------------------------------------------------------

TRIBITS_ADD_LIBRARY(
  kokkosdevicenuma
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS kokkos
  )
