
INCLUDE(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

# TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

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
ReportHandler.hpp
RuntimeDoomed.hpp
RuntimeMessage.hpp
RuntimeWarning.hpp
WallTime.hpp
FormatTime.hpp
FormatMemorySize.hpp
  )

IF ( STK_ENABLE_BoostLib )
  APPEND_SET(HEADERS ProgramOptions.hpp)
ENDIF()

APPEND_SET(SOURCES
CPUTime.cpp
Demangle.cpp
LogControl.cpp
OutputLog.cpp
ProductRegistry.cpp
ReportHandler.cpp
RuntimeDoomed.cpp
RuntimeMessage.cpp
RuntimeWarning.cpp
WallTime.cpp
FormatTime.cpp
FormatMemorySize.cpp
  )

IF ( STK_ENABLE_BoostLib )
  APPEND_SET(SOURCES ProgramOptions.cpp)
ENDIF()


IF (STK_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  stk_util_env
  DEPLIBS stk_util_util stk_util_parallel
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

INSTALL(FILES ${HEADERS} DESTINATION
        ${CMAKE_INSTALL_PREFIX}/${Trilinos_INSTALL_INCLUDE_DIR}/stk_util/environment/)
