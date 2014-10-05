INCLUDE(TribitsLibraryMacros)

SET(HEADERS "")
SET(SOURCES "")

#
# Belos
#

IF (${PACKAGE_NAME}_ENABLE_Belos)

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/belos)

  APPEND_SET(HEADERS
    belos/BelosOperatorT.hpp
    belos/BelosXpetraAdapter.hpp
    belos/BelosXpetraAdapterMultiVector.hpp
    belos/BelosXpetraAdapterOperator.hpp
    belos/BelosMueLuAdapter.hpp
    )

ENDIF()

#
# Epetra
#

IF (${PACKAGE_NAME}_ENABLE_Epetra)

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/epetra)

  APPEND_SET(HEADERS
    epetra/MueLu_EpetraOperator.hpp
    epetra/MueLu_CreateEpetraPreconditioner.hpp
    )

  APPEND_SET(SOURCES
    epetra/MueLu_EpetraOperator.cpp
    )

ENDIF()

#
# Tpetra
#

#IF (${PACKAGE_NAME}_ENABLE_Tpetra)

#  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/tpetra)

#  APPEND_SET(HEADERS
#    tpetra/MueLu_TpetraOperator.hpp
#    )

#  APPEND_SET(SOURCES
#    tpetra/MueLu_TpetraOperator.cpp
#    )

#ENDIF()
IF (${PACKAGE_NAME}_ENABLE_Tpetra)

   SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/tpetra)

   APPEND_GLOB(HEADERS ${DIR}/*.hpp)
   APPEND_GLOB(SOURCES ${DIR}/*.cpp)

   TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

   IF(${PACKAGE_NAME}_ENABLE_Stratimikos AND ${PACKAGE_NAME}_ENABLE_ThyraTpetraAdapters)
     SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/stratimikos-tpetra)
     APPEND_GLOB(HEADERS ${DIR}/*.hpp)
     APPEND_GLOB(SOURCES ${DIR}/*.cpp)
   ENDIF()

ENDIF()

#
# Stratimikos
#

IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Stratimikos)

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/stratimikos/src)

  APPEND_SET(HEADERS
    stratimikos/src/Thyra_MueLuPreconditionerFactory.hpp
    )

  APPEND_SET(SOURCES
    stratimikos/src/Thyra_MueLuPreconditionerFactory.cpp
    )

ENDIF()

#
#
#

# Must glob the binary dir last to get all of the auto-generated headers
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

TRIBITS_ADD_LIBRARY(
  muelu-adapters
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS muelu
  )

# Set linker language explicitly for CUDA builds.
set_target_properties(muelu-adapters PROPERTIES LINKER_LANGUAGE CXX)
