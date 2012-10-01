# CMakeLists inspired by Tpetra/src/CMakeLists.txt

INCLUDE(TribitsLibraryMacros)
INCLUDE(TrilinosCreateClientTemplateHeaders)

# Note: cannot use the macro TRIBITS_CONFIGURE_FILE() because 
# Xpetra_config.hpp.in is not in ${PACKAGE_SOURCE_DIR}/cmake
## TRIBITS_CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Xpetra_config.hpp)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")

SET(SOURCES "mac.cpp")
IF (${PACKAGE_NAME}_ENABLE_Epetra)
  APPEND_GLOB(SOURCES 
    StridedMap/Xpetra_StridedEpetraMap.cpp
    )
ENDIF()

#SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})

IF (NOT Trilinos_ENABLE_Kokkos)
  SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/FakeKokkos)
  APPEND_GLOB(HEADERS FakeKokkos/*.hpp)
ENDIF()

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Matrix)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/StridedMap)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/BlockedCrsMatrix)

#The HEADERS variable is the list of headers that will be copied to the install directory during "make install".
#Note: Do not use */*.hpp. Append FakeKokkos/*.hpp only if Kokkos is disabled.
APPEND_GLOB(HEADERS Matrix/*.hpp)
APPEND_GLOB(HEADERS BlockedCrsMatrix/*.hpp)
APPEND_GLOB(HEADERS StridedMap/*.hpp)

TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

# Must glob the binary dir last to get all of the auto-generated headers
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
#APPEND_SET(HEADERS ../${DIR}/Xpetra_config.hpp)

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  xpetra-sup
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS xpetra
  )
