
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
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/stk_interface)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/evaluators)

# FAD
APPEND_SET(HEADERS
  Panzer_STK_config.hpp
  Panzer_STK_Version.hpp
  Panzer_STK_ConfigDefs.hpp
  Panzer_STKConnManager.hpp
  Panzer_STK_SetupUtilities.hpp
  Panzer_STK_SetupUtilitiesT.hpp
  stk_interface/Panzer_STK_Interface.hpp
  stk_interface/Panzer_STK_MeshFactory.hpp
  stk_interface/Panzer_STK_SquareQuadMeshFactory.hpp
  stk_interface/Panzer_STK_ExodusReaderFactory.hpp
  stk_interface/Panzer_STK_CubeHexMeshFactory.hpp
  stk_interface/Panzer_STK_MultiBlockMeshFactory.hpp
  stk_interface/Panzer_STK_PeriodicBC_Matcher.hpp
  stk_interface/Panzer_STK_PeriodicBC_MatcherT.hpp
  stk_interface/Panzer_STK_PeriodicBC_Parser.hpp
  stk_interface/Panzer_STK_PeriodicBC_MatchConditions.hpp
  evaluators/Panzer_STK_GatherFields.hpp
  evaluators/Panzer_STK_GatherFieldsT.hpp
  Panzer_STK_Utilities.hpp
  ) 

APPEND_SET(SOURCES
  Panzer_STK_Version.cpp
  Panzer_STKConnManager.cpp
  Panzer_STK_SetupUtilities.cpp
  stk_interface/Panzer_STK_Interface.cpp
  stk_interface/Panzer_STK_SquareQuadMeshFactory.cpp
  stk_interface/Panzer_STK_ExodusReaderFactory.cpp
  stk_interface/Panzer_STK_CubeHexMeshFactory.cpp
  stk_interface/Panzer_STK_MultiBlockMeshFactory.cpp
  stk_interface/Panzer_STK_PeriodicBC_Matcher.cpp
  stk_interface/Panzer_STK_PeriodicBC_Parser.cpp
  Panzer_STK_Utilities.cpp
  )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  panzer_stk
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
