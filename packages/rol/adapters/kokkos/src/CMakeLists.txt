

#
# A) Package-specific configuration options
#

# TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# All headers
#

#
# vector headers
#

IF (ROL_ENABLE_Kokkos) 
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/vector)
  APPEND_SET(HEADERS
    vector/ROL_KokkosVector.hpp
    vector/ROL_KokkosVector_decl.hpp
    vector/ROL_KokkosVector_def.hpp
  )
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/elementwise)
  APPEND_SET(HEADERS
   elementwise/ROL_Kokkos_UnaryFunctions_decl.hpp
   elementwise/ROL_Kokkos_UnaryFunctions_def.hpp
  )
#  SET_AND_INC_DIRS( DIR ${CMAKE_CURRENT_SOURCE_DIR}/elementwise )
#  SET_AND_INC_DIRS( DIR ${CMAKE_CURRENT_SOURCE_DIR}/vector )
#  INCLUDE_DIRECTORIES( ${DIR} )
#  APPEND_GLOB(HEADERS ${DIR}/*.hpp)
ENDIF()



INSTALL(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})
