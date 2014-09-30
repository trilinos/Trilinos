INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(TribitsSubPackageMacros)

TRIBITS_SUBPACKAGE(TPL)

SET(TRILINOS_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})

IF( Kokkos_ENABLE_Cuda )

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

  FOREACH( CUB_DIR
           cub
           cub/block
           cub/block/specializations
           cub/device/
           cub/device/block
           cub/device/block/specializations
           cub/grid
           cub/host
           cub/thread
           cub/warp
           cub/warp/specializations )
    SET( CUB_HEADER "" )
    FILE( GLOB CUB_HEADER ${CUB_DIR}/*.cuh )
    INSTALL( FILES ${CUB_HEADER} DESTINATION ${TRILINOS_INCDIR}/${CUB_DIR} )

  ENDFOREACH(CUB_DIR)

ENDIF()

# Need this to set up the ${KokkosTPL_INCLUDE_DIRS} variable
# required by downstream subpackages.
#
# wish: TRIBITS_ADD_HEADER_LIBRARY
#
TRIBITS_ADD_LIBRARY(
    kokkosTPL_unused_dummy
    HEADERS ""
    SOURCES KokkosTPL_dummy.cpp
    DEPLIBS 
    )

TRIBITS_SUBPACKAGE_POSTPROCESS()

