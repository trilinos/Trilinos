#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")

# To include the configured header above
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  SimpleCxx_HelloWorld.hpp
  )

APPEND_SET(SOURCES
  SimpleCxx_HelloWorld.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  simplecxx
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

#
# D) Define an executable that gets installed
#

TRIBITS_ADD_EXECUTABLE(simplecxx-helloworld  NOEXEPREFIX  NOEXESUFFIX
  SOURCES  SimpleCxx_HelloWorld_Prog.cpp
  INSTALLABLE
  )
