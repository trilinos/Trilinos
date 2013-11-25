INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddOptionAndDefine)

TRIBITS_CONFIGURE_FILE( ${PACKAGE_NAME}_config.h )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS "")
SET(SOURCES "")

FILE(GLOB HEADERS *.hpp)
FILE(GLOB SOURCES *.cpp)
# Make sure that the subpackage's automatically generated 
# configuration header file gets installed, along with the
# human-authored header files.  PACKAGE_NAME here indicates
# the subpackage's name.
LIST( APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h )

# Define targets for the subpackage's library.
# No need to use DEPLIBS, since this subpackage will only make one library;
# the subpackage's dependencies handle dependencies outside of the subpackage.
TRIBITS_ADD_LIBRARY(
    kokkoscompat
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    )

IF (Kokkos_ENABLE_Cuda)
  FILE(GLOB CUDA_SOURCES *.cu)
  PACKAGE_ADD_CUDA_LIBRARY(kokkoscompat_cuda SOURCES ${CUDA_SOURCES})
ENDIF()

SET(TRILINOS_INCDIR ${CMAKE_INSTALL_PREFIX}/${Trilinos_INSTALL_INCLUDE_DIR})

INSTALL(FILES ${HEADERS} DESTINATION ${TRILINOS_INCDIR}/)

#-----------------------------------------------------------------------------
