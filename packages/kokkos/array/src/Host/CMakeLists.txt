INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddOptionAndDefine)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS "")
SET(SOURCES "")

FILE(GLOB HEADERS KokkosArray_*.hpp)

#-----------------------------------------------------------------------------
# Sources always included

LIST(APPEND SOURCES KokkosArray_Host_Impl.cpp KokkosArray_HostSpace.cpp)

#-----------------------------------------------------------------------------
# Hardware locality detection functionality depends upon TPL availability

IF( TPL_ENABLE_HWLOC )
#
#  http://www.open-mpi.org/projects/hwloc/
#  Version 1.3 or greater
#
LIST(APPEND SOURCES KokkosArray_Host_hwloc.cpp)
#
ELSE()
LIST(APPEND SOURCES KokkosArray_Host_hwloc_unavailable.cpp)
ENDIF()

#-----------------------------------------------------------------------------
# Thread spawning functionality depend upon TPL availability

IF( TPL_ENABLE_Pthread )
LIST(APPEND SOURCES KokkosArray_Host_pthread.cpp)
ELSEIF( TPL_ENABLE_Winthread )
LIST(APPEND SOURCES KokkosArray_Host_winthread.cpp)
ELSE()
LIST(APPEND SOURCES KokkosArray_Host_nothread.cpp)
ENDIF()

#-----------------------------------------------------------------------------

TRIBITS_ADD_LIBRARY(
  kokkosarray_devicehost
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS 
  )

SET(TRILINOS_INCDIR ${CMAKE_INSTALL_PREFIX}/${Trilinos_INSTALL_INCLUDE_DIR})

INSTALL(FILES ${HEADERS} DESTINATION ${TRILINOS_INCDIR}/Host/)


