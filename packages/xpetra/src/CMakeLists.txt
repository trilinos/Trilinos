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

############################################################################################
FUNCTION(XPETRA_PROCESS_ETI_TEMPLATE ETI_CLASSES TEMPLATE_FILE PROCESSED_FILE SOURCES_LIST)
  SET(SRCS "")
  FOREACH(CLASS ${ETI_CLASSES})
    # find spearator of class name and (optional) conditions
    string(FIND ${CLASS} "-" DASHPOS)
    IF ('${DASHPOS}' STREQUAL '-1')
      string(SUBSTRING ${CLASS} 0 ${DASHPOS} CLASS_NAME)
      #MESSAGE("New class name = ${CLASS_NAME}")
      set (CONDITION_NAME "")
      set (CONDITION_NAME_END "")      
    ELSE()
      string(SUBSTRING ${CLASS} 0 ${DASHPOS} CLASS_NAME)
      #MESSAGE("New class name = ${CLASS_NAME}")
      string(SUBSTRING ${CLASS} ${DASHPOS} -1 CONDITION_NAME)
      string(REPLACE "[" "(" CONDITION_NAME ${CONDITION_NAME})
      string(REPLACE "]" ")" CONDITION_NAME ${CONDITION_NAME})
      string(REPLACE "." " " CONDITION_NAME ${CONDITION_NAME})
      string(REPLACE "-" "" CONDITION_NAME ${CONDITION_NAME})
      string(REPLACE "?" "#" CONDITION_NAME ${CONDITION_NAME})
      string(STRIP CONDITION_NAME ${CONDITION_NAME})
      #MESSAGE("New condition name = ${CONDITION_NAME}")
      set (CONDITION_NAME_END "#endif")
    ENDIF()
    
    string(REPLACE "::" "_" CLASS_FILE_NAME "${CLASS_NAME}")
    string(TOUPPER "${CLASS_FILE_NAME}" UPPER_CASE_CLASS)
    string(REPLACE "CLASS_FILE_NAME" "${CLASS_FILE_NAME}" FINAL_FILE_NAME "${PROCESSED_FILE}")
    
    # the following generates one cpp file for all instantiations and 
    # enabled configurations
    CONFIGURE_FILE(${TEMPLATE_FILE} ${FINAL_FILE_NAME})
    #MESSAGE("Configure file ${FINAL_FILE_NAME} using template ${TEMPLATE_FILE}")
    SET(SRCS ${SRCS} ${FINAL_FILE_NAME})
    
  ENDFOREACH()
  SET(${SOURCES_LIST} ${SRCS} PARENT_SCOPE)
ENDFUNCTION(XPETRA_PROCESS_ETI_TEMPLATE)

##
## XPETRA_ETI_CPP_SOURCES collects all filenames of the Xpetra ETI *.cpp files
## These are later added to the gloabl list of source files (SOURCES)
##
SET(SOURCES "")
SET(XPETRA_ETI_CPP_SOURCES "")

##
## only build ETI *.cpp files if ETI is enabled.
##
IF(${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION)
  ##
  ## XPETRA_SC_LO_GO_NO_ETI_CLASSES is a list of all class names which ETI 
  ## *.cpp files have to be generated for.
  ##
  GLOBAL_SET(XPETRA_SC_LO_GO_NO_ETI_CLASSES  )
  INCLUDE(Utils/ExplicitInstantiation/ETI_SC_LO_GO_NO_classes.cmake)

  XPETRA_PROCESS_ETI_TEMPLATE(
    "${XPETRA_SC_LO_GO_NO_ETI_CLASSES}"
    Utils/ExplicitInstantiation/Xpetra_ETI_SC_LO_GO_NO.tmpl
    "Utils/ExplicitInstantiation/CLASS_FILE_NAME.cpp"
    Xpetra_ETI_SC_LO_GO_NO_SRC)
  LIST(APPEND XPETRA_ETI_CPP_SOURCES ${Xpetra_ETI_SC_LO_GO_NO_SRC})

  MESSAGE("XPETRA_ETI_CPP_SOURCES=${XPETRA_ETI_CPP_SOURCES}")

  # At this point, XPETRA_ETI_CPP_SOURCES contains the Xpetra classes
  # that will be ETI'd.

  APPEND_GLOB(SOURCES ${XPETRA_ETI_CPP_SOURCES})

  SET(SOURCES ${XPETRA_ETI_CPP_SOURCES})
  # grab any hand rolled ETI source files
  APPEND_GLOB(SOURCES ${DIR}/Utils/ExplicitInstantiation/*.cpp)
  
ENDIF()

############################################################################################

SET(HEADERS "")

APPEND_GLOB(SOURCES
  Utils/Xpetra_Utils.cpp
  Utils/Xpetra_EpetraUtils.cpp
  Platform/Xpetra_DefaultPlatform.cpp
  )
IF (${PACKAGE_NAME}_ENABLE_Epetra)
  APPEND_GLOB(SOURCES
    Map/Xpetra_EpetraMap.cpp
    Import/Xpetra_EpetraImport.cpp
    Export/Xpetra_EpetraExport.cpp
    MultiVector/Xpetra_EpetraMultiVector.cpp
    MultiVector/Xpetra_EpetraIntMultiVector.cpp
    CrsGraph/Xpetra_EpetraCrsGraph.cpp
    CrsMatrix/Xpetra_EpetraCrsMatrix.cpp
    Vector/Xpetra_EpetraVector.cpp
    Vector/Xpetra_EpetraIntVector.cpp
    )
ENDIF()

MESSAGE("JHU: CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}")
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/CrsMatrix NOSIERRABJAM)

#SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})

# mfh 18 Dec 2014: The "real" Kokkos(Classic) is now the TpetraClassic
# subpackage of Tpetra.  Thus, we have to test whether Tpetra is
# enabled, not whether Kokkos is enabled.  Furthermore, it wasn't
# previously correct to test if Kokkos was enabled, since Kokkos could
# have been enabled without the KokkosClassic subpackage having been
# enabled.
#
#IF (NOT ${PROJECT_NAME}_ENABLE_Kokkos)
IF (NOT ${PROJECT_NAME}_ENABLE_Tpetra)
  SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/FakeKokkos)
  APPEND_GLOB(HEADERS FakeKokkos/*.hpp)
ENDIF()

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/CrsGraph)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/CrsMatrix)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/BlockedCrsMatrix)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/DistObject)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Export)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Headers)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Import)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Map)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/BlockedMap)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/MultiVector)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/BlockedMultiVector)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/BlockedVector)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Operator)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Platform)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/RowGraph)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/RowMatrix)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Utils)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ForwardDeclaration)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Vector)

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
APPEND_GLOB(HEADERS BlockedMap/*.hpp)
APPEND_GLOB(HEADERS MultiVector/*.hpp)
APPEND_GLOB(HEADERS BlockedMultiVector/*.hpp)
APPEND_GLOB(HEADERS BlockedVector/*.hpp)
APPEND_GLOB(HEADERS Operator/*.hpp)
APPEND_GLOB(HEADERS Platform/*.hpp)
APPEND_GLOB(HEADERS RowGraph/*.hpp)
APPEND_GLOB(HEADERS RowMatrix/*.hpp)
APPEND_GLOB(HEADERS Utils/ForwardDeclaration/*.hpp)
APPEND_GLOB(HEADERS Utils/*.hpp)
APPEND_GLOB(HEADERS Vector/*.hpp)

#TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

# Must glob the binary dir last to get all of the auto-generated headers
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  xpetra
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
