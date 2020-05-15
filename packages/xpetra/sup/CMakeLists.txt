# CMakeLists inspired by Tpetra/src/CMakeLists.txt

INCLUDE(TrilinosCreateClientTemplateHeaders)

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

#
# SRC
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

  # Include the function XPETRA_PROCESS_ETI_TEMPLATE
  include(${CMAKE_CURRENT_SOURCE_DIR}/../src/Utils/ExplicitInstantiation/xpetra-process-eti-template.cmake)


  ##
  ## XPETRA_LO_GO_NO_ETI_CLASSES is a list of all calss names which ETI
  ## *.cpp files are generated for using the 3-arg template (LO-GO-NO)
  ##
  GLOBAL_SET( XPETRA_LO_GO_NO_ETI_CLASSES )

  # Include the list of all Xpetra classes that should be ETI'd 
  INCLUDE(Utils/ExplicitInstantiation/ETI_LO_GO_NO_classes.cmake)

  # generate the cpp files and record the file names
  XPETRA_PROCESS_ETI_TEMPLATE(
    "${XPETRA_LO_GO_NO_ETI_CLASSES}"
    ${Trilinos_SOURCE_DIR}/packages/xpetra/src/Utils/ExplicitInstantiation/Xpetra_ETI_LO_GO_NO.tmpl
    "Utils/ExplicitInstantiation/CLASS_FILE_NAME.cpp"
    Xpetra_ETI_LO_GO_NO_SRC)

  LIST(APPEND XPETRA_ETI_CPP_SOURCES ${Xpetra_ETI_LO_GO_NO_SRC})


  ##
  ## XPETRA_SC_LO_GO_NO_ETI_CLASSES is a list of all class names which ETI 
  ## *.cpp files have to be generated for using the 4-arg template (SC-LO-GO-NO)
  ##
  GLOBAL_SET(XPETRA_SC_LO_GO_NO_ETI_CLASSES)

  # include the list of all Xpetra classes that should be ETI'd
  INCLUDE(Utils/ExplicitInstantiation/ETI_SC_LO_GO_NO_classes.cmake)

  # generate the cpp files and record the file names
  XPETRA_PROCESS_ETI_TEMPLATE(
    "${XPETRA_SC_LO_GO_NO_ETI_CLASSES}"
    ${Trilinos_SOURCE_DIR}/packages/xpetra/src/Utils/ExplicitInstantiation/Xpetra_ETI_SC_LO_GO_NO.tmpl
    "Utils/ExplicitInstantiation/CLASS_FILE_NAME.cpp"
    Xpetra_ETI_SC_LO_GO_NO_SRC)
  LIST(APPEND XPETRA_ETI_CPP_SOURCES ${Xpetra_ETI_SC_LO_GO_NO_SRC})


  # At this point, XPETRA_ETI_CPP_SOURCES contains the Xpetra classes
  # that will be ETI'd.

  APPEND_GLOB(SOURCES ${XPETRA_ETI_CPP_SOURCES})

  SET(SOURCES ${XPETRA_ETI_CPP_SOURCES})
  # grab any hand rolled ETI source files
  APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ExplicitInstantiation/*.cpp)
  
ENDIF() # if XPETRA_ENABLE_EXPLICIT_INSTANTIATION


SET(HEADERS "")
APPEND_GLOB(HEADERS */*.hpp)
APPEND_GLOB(SOURCES */*.cpp)

# FIXME JHU is this still needed?
APPEND_GLOB(SOURCES mac.cpp) # Avoid empty library (MacOs 'ar' utility doesn't work on empty library).

# Create the .hpp file in BUILD/packages/xpetra/sup
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/Matrix     NOSIERRABJAM)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/StridedMap NOSIERRABJAM)

#
#
#

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Matrix)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/StridedMap)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Utils)

TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

# Must glob the binary dir last to get all of the auto-generated headers
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#
# Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  xpetra-sup
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS xpetra
  )
