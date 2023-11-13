

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

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# All headers
#

#
# sol headers
#

IF (ROL_ENABLE_TriKota)

  TRIBITS_SET_AND_INC_DIRS( DIR ${CMAKE_CURRENT_SOURCE_DIR}/sol )
  TRIBITS_INCLUDE_DIRECTORIES( ${DIR} )
  APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#  TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/sol)
#  TRIBITS_INCLUDE_DIRECTORIES(${DIR}
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

INSTALL(FILES ${HEADERS} DESTINATION ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})
