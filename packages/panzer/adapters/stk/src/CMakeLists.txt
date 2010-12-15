
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.hpp
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/albany)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/stk_interface)

# FAD
APPEND_SET(HEADERS
  Panzer_STK_config.hpp
  Panzer_STK_Version.hpp
  Panzer_STK_ConfigDefs.hpp
#  Panzer_AlbanySTKConnManager.hpp
#  albany/Albany_AbstractSTKMeshStruct.hpp
#  albany/Albany_Quad2DSTKMeshStruct.hpp
#  albany/Albany_AbstractDiscretization.hpp
#  albany/Albany_STKDiscretization.hpp
  Panzer_STKConnManager.hpp
  Panzer_STK_SetupUtilities.hpp
  Panzer_STK_SetupUtilitiesT.hpp
  stk_interface/Panzer_STK_Interface.hpp
  stk_interface/Panzer_STK_MeshFactory.hpp
  stk_interface/Panzer_STK_SquareQuadMeshFactory.hpp
  stk_interface/Panzer_STK_ExodusReaderFactory.hpp
  ) 

APPEND_SET(SOURCES
  Panzer_STK_Version.cpp
  #Panzer_AlbanySTKConnManager.cpp
  #albany/Albany_Quad2DSTKMeshStruct.cpp
  #albany/Albany_STKDiscretization.cpp
  Panzer_STKConnManager.cpp
  Panzer_STK_SetupUtilities.cpp
  stk_interface/Panzer_STK_Interface.cpp
  stk_interface/Panzer_STK_SquareQuadMeshFactory.cpp
  stk_interface/Panzer_STK_ExodusReaderFactory.cpp
  )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  panzer_stk
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
