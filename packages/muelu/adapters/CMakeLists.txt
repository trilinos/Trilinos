INCLUDE(TribitsLibraryMacros)

SET(HEADERS "")
SET(SOURCES "")

#
# Belos
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/belos)

IF (${PACKAGE_NAME}_ENABLE_Belos)

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
    )

  APPEND_SET(SOURCES
    epetra/MueLu_EpetraOperator.cpp
    )

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

APPEND_SET(SOURCES mac.cpp)

TRIBITS_ADD_LIBRARY(
  muelu-adapters
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS muelu
  )
