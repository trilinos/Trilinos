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
IF (${PACKAGE_NAME}_ENABLE_Tpetra)

   SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/tpetra)

   APPEND_GLOB(HEADERS ${DIR}/*.hpp)
   APPEND_GLOB(SOURCES ${DIR}/*.cpp)

   TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

   IF (${PACKAGE_NAME}_ENABLE_Stratimikos AND ${PACKAGE_NAME}_ENABLE_ThyraTpetraAdapters)
     SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/stratimikos-tpetra)
     APPEND_SET(SOURCES
       ${DIR}/Stratimikos_MueluTpetraHelpers.cpp
       )
     APPEND_SET(HEADERS
       ${DIR}/Stratimikos_MueluTpetraHelpers.hpp
       ${DIR}/Thyra_MueLuTpetraPreconditionerFactory_decl.hpp
       ${DIR}/Thyra_MueLuTpetraPreconditionerFactory_def.hpp
       ${DIR}/Thyra_MueLuTpetraPreconditionerFactory.hpp
       )
     IF (${PACKAGE_NAME}_ENABLE_Teko)
       APPEND_SET(HEADERS
         ${DIR}/Thyra_MueLuTpetraQ2Q1PreconditionerFactory_decl.hpp
         ${DIR}/Thyra_MueLuTpetraQ2Q1PreconditionerFactory_def.hpp
         ${DIR}/Thyra_MueLuTpetraQ2Q1PreconditionerFactory.hpp
         )
     ENDIF()
   ENDIF()

ENDIF()

#
# Stratimikos
#
IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Stratimikos)

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/stratimikos)

  APPEND_SET(HEADERS
    stratimikos/Thyra_MueLuPreconditionerFactory.hpp
    )

  APPEND_SET(SOURCES
    stratimikos/Thyra_MueLuPreconditionerFactory.cpp
    )

ENDIF()

# Must glob the binary dir last to get all of the auto-generated headers
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

TRIBITS_ADD_LIBRARY(
  muelu-adapters
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS muelu
  ADDED_LIB_TARGET_NAME_OUT MUELU_ADAPTERS_LIBNAME
  )

# Set linker language explicitly for CUDA builds.
set_target_properties(
  ${MUELU_ADAPTERS_LIBNAME}
  PROPERTIES LINKER_LANGUAGE CXX
  )
