
INCLUDE(PackageLibraryMacros)

#
# Define the configure and build of Dakota libraries
#

# ToDo: Uncomment this below

#EXECUTE_PROCESS(
#  WORKING_DIRECTORY ${PACKAGE_SOURCE_DIR}/Dakota
#  COMMAND ./configure
#    --without-graphics CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
#    --prefix=${PACKAGE_BINARY_DIR}/Dakota
#    --with-teuchos-include=${PACKAGE_SOURCE_DIR}/../teuchos/src
#    --with-teuchos-lib=${PACKAGE_BINARY_DIR}/../teuchos/src
#    --without-graphics --without-acro --without-ddace   --without-evidence
#    --without-fsudace --without-appspack  --without-jega --without-nlpql
#    --without-optpp --without-psuade --disable-f90 --disable-mpi
#  OUTPUT_FILE ${PACKAGE_SOURCE_DIR}/Dakota/configure.out
#  )
#
## ToDo: Add *all* compiler options.  See Makefile.export.in and
## PackageWriteExportMakefile.cmake for how to grab *all* the compiler options
## that you need to be consistent with building against Trilinos.
#
#ADD_CUSTOM_COMMAND(
#  OUTPUT ${PACKAGE_SOURCE_DIR}/Dakota/src/dakota
#  COMMAND make # ToDo: Get right varible
#  WORKING_DIRECTORY ${PACKAGE_SOURCE_DIR}/Dakota
#  )
#
#ADD_CUSTOM_TARGET( ${PACKAGE_NAME}_Dakota_libs ALL
#  DEPENDS ${PACKAGE_SOURCE_DIR}/Dakota/src/dakota
#  )

#
# Build the TriKota library
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# FAD
APPEND_SET(HEADERS
    TriKota_ConfigDefs.hpp
    TriKota_DirectApplicInterface.hpp
    TriKota_Driver.hpp
    TriKota_Version.hpp
  ) 

APPEND_SET(SOURCES
    TriKota_DirectApplicInterface.cpp
    TriKota_Driver.cpp
    TriKota_Version.cpp
  )

SET(DAKOTA_LIBS -ldakota -lblob)

# ToDo: Above: Set up DAKOTA_LIBS correctly!

PACKAGE_ADD_LIBRARY(
  trikota
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  #DEPLIBS ${DAKOTA_LIBS}
  )

# ToDo: Above: Uncomment DEPLIBS
