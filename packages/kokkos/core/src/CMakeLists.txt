INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddOptionAndDefine)

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#-----------------------------------------------------------------------------

SET(HEADERS "")
SET(SOURCES "")

FILE(GLOB HEADERS Kokkos*.hpp)

SET(TRILINOS_INCDIR ${CMAKE_INSTALL_PREFIX}/${Trilinos_INSTALL_INCLUDE_DIR})

INSTALL(FILES ${HEADERS} DESTINATION ${TRILINOS_INCDIR}/)

#-----------------------------------------------------------------------------

ADD_SUBDIRECTORY(Host)
ADD_SUBDIRECTORY(Threads)
ADD_SUBDIRECTORY(impl)

IF(Kokkos_ENABLE_OpenMP)
  ADD_SUBDIRECTORY(OpenMP)
ENDIF()

IF(Kokkos_ENABLE_Cuda)
  ADD_SUBDIRECTORY(Cuda)
ENDIF()

