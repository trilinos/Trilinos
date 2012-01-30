
INCLUDE(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

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
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/evaluators)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/responses)

# auto generate headers, add to dependencies and include directories
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR})
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/dofmngr)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/evaluators)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/responses)
APPEND_GLOB(HEADERS ${CMAKE_CURRENT_BINARY_DIR}/*.hpp)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

APPEND_GLOB(HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
APPEND_GLOB(HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/evaluators/*.hpp)
APPEND_GLOB(HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/responses/*.hpp)

APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/evaluators/*.cpp)
APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/responses/*.cpp)

APPEND_SET(HEADERS
  dofmngr/Panzer_FieldPattern.hpp
  dofmngr/Panzer_IntrepidFieldPattern.hpp
  dofmngr/Panzer_GeometricAggFieldPattern.hpp
  dofmngr/Panzer_FieldAggPattern.hpp
  dofmngr/Panzer_ConnManager.hpp
  dofmngr/Panzer_DOFManager.hpp
  dofmngr/Panzer_DOFManagerT.hpp
  dofmngr/Panzer_BlockedDOFManager.hpp
  dofmngr/Panzer_BlockedDOFManagerT.hpp
  dofmngr/Panzer_UniqueGlobalIndexer_Utilities.hpp
  dofmngr/Panzer_UniqueGlobalIndexer_UtilitiesT.hpp
  dofmngr/Panzer_UniqueGlobalIndexer.hpp
  dofmngr/Panzer_EpetraLinearObjFactory.hpp
  dofmngr/Panzer_EpetraLinearObjFactoryT.hpp
  dofmngr/Panzer_EpetraLinearObjContainer.hpp
  dofmngr/Panzer_UniqueGlobalIndexerFactory.hpp
  dofmngr/Panzer_DOFManagerFactory.hpp
  dofmngr/Panzer_DOFManagerFactoryT.hpp
  dofmngr/Panzer_LinearObjFactory.hpp
  ) 

APPEND_SET(SOURCES
  dofmngr/Panzer_FieldPattern.cpp
  dofmngr/Panzer_IntrepidFieldPattern.cpp
  dofmngr/Panzer_GeometricAggFieldPattern.cpp
  dofmngr/Panzer_FieldAggPattern.cpp
  dofmngr/Panzer_DOFManager.cpp
  dofmngr/Panzer_BlockedDOFManager.cpp
  dofmngr/Panzer_EpetraLinearObjFactory.cpp
  )

IF(${PACKAGE_NAME}_ENABLE_Stokhos)
  APPEND_SET(HEADERS
    dofmngr/Panzer_SGEpetraLinearObjFactory.hpp
    dofmngr/Panzer_SGEpetraLinearObjFactoryT.hpp
    dofmngr/Panzer_SGEpetraLinearObjContainer.hpp
    )
  APPEND_SET(SOURCES
    dofmngr/Panzer_SGEpetraLinearObjFactory.cpp
    dofmngr/Panzer_SGEpetraLinearObjContainer.cpp
    )
ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  panzer
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

