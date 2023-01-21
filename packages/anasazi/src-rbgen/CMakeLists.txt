

#
# A) Package-specific configuration options
#

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(NOINSTALLHEADERS "")
SET(SOURCES "")

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# Abstract headers
APPEND_SET(HEADERS
  RBGen_ConfigDefs.h
  RBGen_FileIOFactory.hpp
  RBGen_FileIOHandler.hpp
  RBGen_Filter.hpp
  RBGen_MethodFactory.hpp
  RBGen_Method.hpp
  RBGen_NoPreprocessor.hpp
  RBGen_Params.h
  RBGen_PODMethod.hpp
  RBGen_PreprocessorFactory.hpp
  RBGen_Preprocessor.hpp
  RBGen_Version.h
  )

APPEND_SET(SOURCES
  RBGen_Params.cpp
  )

#
# C) Define the targets for package's library(s)
#

IF (${PACKAGE_NAME}_ENABLE_RBGen)
  TRIBITS_ADD_LIBRARY(
    rbgen
    HEADERS ${HEADERS}
    NOINSTALLHEADERS ${NOINSTALLHEADERS}
    SOURCES ${SOURCES}
    DEPLIBS anasazi
    )
ENDIF()
