# CMakeLists inspired by Tpetra/src/CMakeLists.txt

INCLUDE(TrilinosCreateClientTemplateHeaders)

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

#
# B) Define the header and source files (and directories)
#

#
# src
#

##
## XPETRA_ETI_CPP_SOURCES collects all filenames of the Xpetra ETI *.cpp files
## These are later added to the global list of source files (SOURCES)
##
SET(SOURCES "")
SET(XPETRA_ETI_CPP_SOURCES "")

##
## only build ETI *.cpp files if ETI is enabled.
##
IF(${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION)
  ##
  ## XPETRA_(Something)_ETI_CLASSES is a list of all class names which ETI
  ## *.cpp files have to be generated for.
  ##

  # include the function XPETRA_PROCESS_ETI_TEMPLATE
  include(Utils/ExplicitInstantiation/xpetra-process-eti-template.cmake)

  GLOBAL_SET(XPETRA_LO_GO_NO_ETI_CLASSES  )
  INCLUDE(Utils/ExplicitInstantiation/ETI_LO_GO_NO_classes.cmake)
  XPETRA_PROCESS_ETI_TEMPLATE(
    "${XPETRA_LO_GO_NO_ETI_CLASSES}"
    ${Trilinos_SOURCE_DIR}/packages/xpetra/src/Utils/ExplicitInstantiation/Xpetra_ETI_LO_GO_NO.tmpl
    "Utils/ExplicitInstantiation/CLASS_FILE_NAME.cpp"
    Xpetra_ETI_LO_GO_NO_SRC)
  LIST(APPEND XPETRA_ETI_CPP_SOURCES ${Xpetra_ETI_LO_GO_NO_SRC})


  GLOBAL_SET(XPETRA_SC_LO_GO_NO_ETI_CLASSES  )
  INCLUDE(Utils/ExplicitInstantiation/ETI_SC_LO_GO_NO_classes.cmake)
  XPETRA_PROCESS_ETI_TEMPLATE(
    "${XPETRA_SC_LO_GO_NO_ETI_CLASSES}"
    ${Trilinos_SOURCE_DIR}/packages/xpetra/src/Utils/ExplicitInstantiation/Xpetra_ETI_SC_LO_GO_NO.tmpl
    "Utils/ExplicitInstantiation/CLASS_FILE_NAME.cpp"
    Xpetra_ETI_SC_LO_GO_NO_SRC)
  LIST(APPEND XPETRA_ETI_CPP_SOURCES ${Xpetra_ETI_SC_LO_GO_NO_SRC})


  GLOBAL_SET(XPETRA_MULTIVECTOR_ETI_CLASSES  )
  INCLUDE(Utils/ExplicitInstantiation/ETI_MULTIVECTOR_classes.cmake)
  XPETRA_PROCESS_ETI_TEMPLATE(
    "${XPETRA_MULTIVECTOR_ETI_CLASSES}"
    ${Trilinos_SOURCE_DIR}/packages/xpetra/src/Utils/ExplicitInstantiation/Xpetra_ETI_MULTIVECTOR.tmpl
    "Utils/ExplicitInstantiation/CLASS_FILE_NAME.cpp"
    Xpetra_ETI_MULTIVECTOR_SRC)
  LIST(APPEND XPETRA_ETI_CPP_SOURCES ${Xpetra_ETI_MULTIVECTOR_SRC})
  MESSAGE(STATUS "Xpetra: Enabling ETI")
#  MESSAGE("XPETRA_ETI_CPP_SOURCES=${XPETRA_ETI_CPP_SOURCES}")

  # At this point, XPETRA_ETI_CPP_SOURCES contains the Xpetra classes
  # that will be ETI'd.

  APPEND_GLOB(SOURCES ${XPETRA_ETI_CPP_SOURCES})

  SET(SOURCES ${XPETRA_ETI_CPP_SOURCES})
  # grab any hand rolled ETI source files
  APPEND_GLOB(SOURCES ${DIR}/Utils/ExplicitInstantiation/*.cpp)
ELSE()
  MESSAGE(STATUS "Xpetra: Disabling ETI")
ENDIF()

############################################################################################

SET(HEADERS "")

APPEND_GLOB(SOURCES
  Utils/Xpetra_Utils.cpp
  Utils/Xpetra_EpetraUtils.cpp
  Utils/Xpetra_IteratorOps.cpp
  Platform/Xpetra_DefaultPlatform.cpp
  )
IF (${PACKAGE_NAME}_ENABLE_Epetra)
  APPEND_GLOB(SOURCES
    Map/Xpetra_EpetraMap.cpp
    Map/Xpetra_EpetraMapFactory.cpp
    Import/Xpetra_EpetraImport.cpp
    Export/Xpetra_EpetraExport.cpp
    MultiVector/Xpetra_EpetraMultiVector.cpp
    MultiVector/Xpetra_EpetraIntMultiVector.cpp
    MultiVector/Xpetra_EpetraMultiVectorFactory.cpp
    CrsGraph/Xpetra_EpetraCrsGraph.cpp
    CrsMatrix/Xpetra_EpetraCrsMatrix.cpp
    Vector/Xpetra_EpetraVector.cpp
    Vector/Xpetra_EpetraIntVector.cpp
    Vector/Xpetra_EpetraVectorFactory.cpp
    )
ENDIF()

# Create the .hpp file in BUILD/packages/xpetra/src
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/BlockedCrsMatrix   NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/BlockedMap         NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/BlockedMultiVector NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/BlockedVector      NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/CrsGraph           NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/CrsMatrix          NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/Export             NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/Map                NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/Matrix             NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/MultiVector        NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/Utils              NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/StridedMap         NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/Vector             NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/Import             NOSIERRABJAM)


#TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})

# mfh 18 Dec 2014: The "real" Kokkos(Classic) is now the TpetraClassic
# subpackage of Tpetra.  Thus, we have to test whether Tpetra is
# enabled, not whether Kokkos is enabled.  Furthermore, it wasn't
# previously correct to test if Kokkos was enabled, since Kokkos could
# have been enabled without the KokkosClassic subpackage having been
# enabled.
#
#IF (NOT ${PROJECT_NAME}_ENABLE_Kokkos)
IF (NOT ${PROJECT_NAME}_ENABLE_Tpetra)
  TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/FakeKokkos)
  APPEND_GLOB(HEADERS FakeKokkos/*.hpp)
ENDIF()

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/CrsGraph)
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/CrsMatrix)
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/BlockedCrsMatrix)
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/DistObject)
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Export)
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Headers)
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Import)
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Map)
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/BlockedMap)
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/MultiVector)
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/BlockedMultiVector)
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/BlockedVector)
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Matrix)
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Operator)
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Platform)
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/RowGraph)
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/RowMatrix)
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/StridedMap)
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Utils)
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ForwardDeclaration)
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Vector)

APPEND_GLOB(SOURCES BlockedCrsMatrix/*.cpp)

#The HEADERS variable is the list of headers that will be copied to the install directory during "make install".
#Note: Do not use */*.hpp. Append FakeKokkos/*.hpp only if Kokkos is disabled.
APPEND_GLOB(HEADERS CrsGraph/*.hpp)
APPEND_GLOB(HEADERS CrsMatrix/*.hpp)
APPEND_GLOB(HEADERS BlockedCrsMatrix/*.hpp)
APPEND_GLOB(HEADERS DistObject/*.hpp)
APPEND_GLOB(HEADERS Export/*.hpp)
APPEND_GLOB(HEADERS Headers/*.hpp)
APPEND_GLOB(HEADERS Import/*.hpp)
APPEND_GLOB(HEADERS Map/*.hpp)
APPEND_GLOB(HEADERS Matrix/*.hpp)
APPEND_GLOB(HEADERS BlockedMap/*.hpp)
APPEND_GLOB(HEADERS MultiVector/*.hpp)
APPEND_GLOB(HEADERS BlockedMultiVector/*.hpp)
APPEND_GLOB(HEADERS BlockedVector/*.hpp)
APPEND_GLOB(HEADERS Operator/*.hpp)
APPEND_GLOB(HEADERS Platform/*.hpp)
APPEND_GLOB(HEADERS RowGraph/*.hpp)
APPEND_GLOB(HEADERS RowMatrix/*.hpp)
APPEND_GLOB(HEADERS StridedMap/*.hpp)
APPEND_GLOB(HEADERS Utils/ForwardDeclaration/*.hpp)
APPEND_GLOB(HEADERS Utils/*.hpp)
APPEND_GLOB(HEADERS Vector/*.hpp)

#TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

# Must glob the binary dir last to get all of the auto-generated headers
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  xpetra
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
# touch CMakeLists.txt because a new file was created in Utils/ExplicitInstantiation of Utils/ForwardDeclaration
# touch CMakeLists.txt because a new file was created in Utils/ExplicitInstantiation of Utils/ForwardDeclaration
