INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddOptionAndDefine)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS "")
SET(SOURCES "")

FILE(GLOB HEADERS *.hpp)
FILE(GLOB SOURCES *.cpp)

#-----------------------------------------------------------------------------
# Sources always included

TRIBITS_ADD_LIBRARY(
  kokkosarray_deviceopenmp
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS 
  )

SET(TRILINOS_INCDIR ${CMAKE_INSTALL_PREFIX}/${Trilinos_INSTALL_INCLUDE_DIR})

INSTALL(FILES ${HEADERS} DESTINATION ${TRILINOS_INCDIR}/OpenMP/)


