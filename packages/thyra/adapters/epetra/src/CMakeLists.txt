

SET(HEADERS "")
SET(SOURCES "")

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

# This allows to include the ThyraEpetraAdapters_config.h file
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

IF (${ThyraEpetraAdapters_ENABLE_Refactor})
  SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/new_interface)
  SET (HEADERS
    ${DIR}/Thyra_EpetraLinearOp.hpp
    ${DIR}/Thyra_EpetraMultiVector.hpp
    ${DIR}/Thyra_EpetraOperatorWrapper.hpp
    ${DIR}/Thyra_EpetraThyraWrappers.hpp
    ${DIR}/Thyra_EpetraVector.hpp
    ${DIR}/Thyra_EpetraVectorSpace.hpp
  )
  SET (SOURCES
    ${DIR}/Thyra_EpetraLinearOp.cpp
    ${DIR}/Thyra_EpetraMultiVector.cpp
    ${DIR}/Thyra_EpetraOperatorWrapper.cpp
    ${DIR}/Thyra_EpetraThyraWrappers.cpp
    ${DIR}/Thyra_EpetraVector.cpp
    ${DIR}/Thyra_EpetraVectorSpace.cpp
  )
ELSE ()
  SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/old_interface)
  APPEND_GLOB(HEADERS ${DIR}/*.hpp)
  APPEND_GLOB(SOURCES ${DIR}/*.cpp)
ENDIF()

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET (SOURCES
  ${SOURCES}
  ${DIR}/Thyra_EpetraOperatorViewExtractorStd.cpp
  ${DIR}/Thyra_get_Epetra_Operator.cpp
)

SET (HEADERS
  ${HEADERS}
  ${DIR}/Thyra_EpetraOperatorViewExtractorBase.hpp
  ${DIR}/Thyra_EpetraOperatorViewExtractorStd.hpp
  ${DIR}/Thyra_EpetraTypes.hpp
  ${DIR}/Thyra_get_Epetra_Operator.hpp
)

TRIBITS_ADD_LIBRARY(
  thyraepetra
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
