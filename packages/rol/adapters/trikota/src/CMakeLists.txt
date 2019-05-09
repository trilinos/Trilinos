

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
# sol headers
#

IF (ROL_ENABLE_TriKota)

  SET_AND_INC_DIRS( DIR ${CMAKE_CURRENT_SOURCE_DIR}/sol )
  INCLUDE_DIRECTORIES( ${DIR} )
  APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#  SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/sol)
#  INCLUDE_DIRECTORIES(${DIR}
#                      ${DIR}/sparse_grids
#  )
#  APPEND_SET(HEADERS ${DIR}/ROL_SparseGridGenerator.hpp
#                     ${DIR}/ROL_SparseGridGeneratorDef.hpp
#                     ${DIR}/sparse_grids/ROL_Quadrature.hpp
#                     ${DIR}/sparse_grids/ROL_QuadratureDef.hpp
#                     ${DIR}/sparse_grids/ROL_Quadrature1dConstructors.hpp
#                     ${DIR}/sparse_grids/ROL_QuadratureTPConstructors.hpp
#                     ${DIR}/sparse_grids/ROL_QuadratureHelpers.hpp
#                     ${DIR}/sparse_grids/ROL_QuadratureTypes.hpp
#  )
ENDIF()

INSTALL(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})
