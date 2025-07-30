TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Graph) # we need this line because this file is not on the header file List of MueLu
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)

# first ensure that these cmake boolean variables are defined

IF(${PACKAGE_NAME}_ENABLE_Epetra
   AND ${PACKAGE_NAME}_ENABLE_EpetraExt
   AND ${PACKAGE_NAME}_ENABLE_Ifpack
   AND ${PACKAGE_NAME}_ENABLE_Amesos
)


ENDIF()
