#
# A) Package-specific configuration options
#

tribits_configure_file(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and include directories)
#

set(HEADERS "")
set(SOURCES "")

tribits_include_directories(${CMAKE_CURRENT_BINARY_DIR})
append_set(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

tribits_include_directories(${CMAKE_CURRENT_SOURCE_DIR})
append_set(HEADERS
  SimpleCxx_HelloWorld.hpp
  )
append_set(SOURCES
  SimpleCxx_HelloWorld.cpp
  )

#
# C) Define the targets for package's library(s)
#

tribits_add_library(
  simplecxx
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

#
# D) Define an executable that gets installed
#

tribits_add_executable(simplecxx-helloworld  NOEXEPREFIX  NOEXESUFFIX
  SOURCES  SimpleCxx_HelloWorld_Prog.cpp
  INSTALLABLE
  )
