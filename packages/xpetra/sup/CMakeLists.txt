# CMakeLists inspired by Tpetra/src/CMakeLists.txt

INCLUDE(TribitsLibraryMacros)
INCLUDE(TrilinosCreateClientTemplateHeaders)

#
# SRC
#

SET(HEADERS "")
APPEND_GLOB(HEADERS */*.hpp)

SET(SOURCES "mac.cpp") # Avoid empty library (MacOs 'ar' utility doesn't work on empty library).
IF (${PACKAGE_NAME}_ENABLE_Epetra)
  APPEND_GLOB(SOURCES 
    StridedMap/Xpetra_StridedEpetraMap.cpp
    )
ENDIF()

#
#
#

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Matrix)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/StridedMap)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/BlockedCrsMatrix)
# Must glob the binary dir last to get all of the auto-generated headers
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})

# unused TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

#
# Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  xpetra-sup
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS xpetra
  )
