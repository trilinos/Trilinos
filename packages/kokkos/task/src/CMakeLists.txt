INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddOptionAndDefine)

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#-----------------------------------------------------------------------------

SET(TRILINOS_INCDIR ${CMAKE_INSTALL_PREFIX}/${Trilinos_INSTALL_INCLUDE_DIR})

#-----------------------------------------------------------------------------

SET(HEADERS_PUBLIC "")
SET(HEADERS_PRIVATE "")
SET(SOURCES "")

FILE(GLOB HEADERS_PUBLIC Kokkos*.hpp)

#-----------------------------------------------------------------------------

FILE(GLOB HEADERS_SERIAL Serial/*.hpp)
FILE(GLOB SOURCES_SERIAL Serial/*.cpp)

LIST(APPEND HEADERS_PRIVATE ${HEADERS_SERIAL} )
LIST(APPEND SOURCES         ${SOURCES_SERIAL} )

INSTALL(FILES ${HEADERS_IMPL} DESTINATION ${TRILINOS_INCDIR}/Serial/)

#-----------------------------------------------------------------------------

# FILE(GLOB HEADERS_THREADS Threads/*.hpp)
# FILE(GLOB SOURCES_THREADS Threads/*.cpp)
# 
# LIST(APPEND HEADERS_PRIVATE ${HEADERS_THREADS} )
# LIST(APPEND SOURCES         ${SOURCES_THREADS} )
# 
# INSTALL(FILES ${HEADERS_THREADS} DESTINATION ${TRILINOS_INCDIR}/Threads/)

#-----------------------------------------------------------------------------

# IF(Kokkos_ENABLE_OpenMP)
# 
#   FILE(GLOB HEADERS_OPENMP OpenMP/*.hpp)
#   FILE(GLOB SOURCES_OPENMP OpenMP/*.cpp)
# 
#   LIST(APPEND HEADERS_PRIVATE ${HEADERS_OPENMP} )
#   LIST(APPEND SOURCES         ${SOURCES_OPENMP} )
# 
#   INSTALL(FILES ${HEADERS_OPENMP} DESTINATION ${TRILINOS_INCDIR}/OpenMP/)
# 
# ENDIF()

#-----------------------------------------------------------------------------

# IF(Kokkos_ENABLE_Cuda)
# 
#  # set nvcc behavior
#
#   SET(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF CACHE BOOL "" FORCE )
#   SET(CUDA_PROPAGATE_HOST_FLAGS              OFF CACHE BOOL "" FORCE )
# 
#   FILE(GLOB HEADERS_CUDA Cuda/*.hpp)
#   FILE(GLOB SOURCES_CUDA Cuda/*.cu)
# 
#   LIST(APPEND HEADERS_PRIVATE ${HEADERS_CUDA} )
# 
#   # Add subpackage include path dependencies for CUDA_COMPILE().
#   CUDA_INCLUDE_DIRECTORIES( ${KokkosTPL_INCLUDE_DIRS} )
#   CUDA_COMPILE(OBJECTS_CUDA ${SOURCES_CUDA})
# 
#   LIST(APPEND SOURCES ${OBJECTS_CUDA} )
# 
#   INSTALL(FILES ${HEADERS_CUDA} DESTINATION ${TRILINOS_INCDIR}/Cuda/)
# 
# ENDIF()

#-----------------------------------------------------------------------------

TRIBITS_ADD_LIBRARY(
    kokkostask
    HEADERS ${HEADERS_PUBLIC}
    NOINSTALLHEADERS ${HEADERS_PRIVATE}
    SOURCES ${SOURCES}
    DEPLIBS
    )

