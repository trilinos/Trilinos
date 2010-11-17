
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
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/dofmngr)

APPEND_SET(HEADERS
  Panzer_ArrayTraits.hpp
  Panzer_Basis.hpp
  Panzer_BasisValues.hpp
  Panzer_BasisValuesT.hpp
  Panzer_BC.hpp
  Panzer_CellData.hpp
  Panzer_ConfigDefs.hpp
  Panzer_Dimension.hpp
  Panzer_InputPhysicsBlock.hpp
  Panzer_IntegrationRule.hpp
  Panzer_IntegrationValues.hpp
  Panzer_IntegrationValuesT.hpp
  Panzer_IntrepidBasisFactory.hpp
  Panzer_Version.hpp
  Panzer_Workset.hpp
  Panzer_Workset_Builder.hpp
  Panzer_Workset_BuilderT.hpp
  dofmngr/Panzer_FieldPattern.hpp
  dofmngr/Panzer_IntrepidFieldPattern.hpp
  dofmngr/Panzer_GeometricAggFieldPattern.hpp
  dofmngr/Panzer_FieldAggPattern.hpp
  dofmngr/Panzer_ConnManager.hpp
  dofmngr/Panzer_DOFManager.hpp
  ) 

APPEND_SET(SOURCES
  Panzer_Basis.cpp
  Panzer_BC.cpp
  Panzer_Dimension.cpp
  Panzer_InputPhysicsBlock.cpp
  Panzer_IntegrationRule.cpp
  Panzer_Version.cpp
  dofmngr/Panzer_FieldPattern.cpp
  dofmngr/Panzer_IntrepidFieldPattern.cpp
  dofmngr/Panzer_GeometricAggFieldPattern.cpp
  dofmngr/Panzer_FieldAggPattern.cpp
  dofmngr/Panzer_DOFManager.cpp
  )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  panzer
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

