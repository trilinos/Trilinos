
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

# PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/..)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../)

#
# Core headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
CPUTime.hpp
Demangle.hpp
LogControl.hpp
OutputLog.hpp
ProductRegistry.hpp
ProgramOptions.hpp
ReportHandler.hpp
RuntimeDoomed.hpp
RuntimeMessage.hpp
RuntimeWarning.hpp
WallTime.hpp
  )

APPEND_SET(SOURCES
CPUTime.cpp
Demangle.cpp
LogControl.cpp
OutputLog.cpp
ProductRegistry.cpp
ProgramOptions.cpp
ReportHandler.cpp
RuntimeDoomed.cpp
RuntimeMessage.cpp
RuntimeWarning.cpp
WallTime.cpp
  )

IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  stk_util_env
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

INSTALL(FILES ${HEADERS} DESTINATION
        ${CMAKE_INSTALL_PREFIX}/include/stk_util/environment/)
