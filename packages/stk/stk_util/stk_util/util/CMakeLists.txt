
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

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS}/../)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS}/detail/)

#
# Core headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
Bootstrap.hpp
ci_string.hpp
ci_traits.hpp
CSet.hpp
FormatTime.hpp
Identifier.hpp
IndentStreambuf.hpp
IndexList.hpp
Marshal.hpp
PairIter.hpp
PrintTable.hpp
SameType.hpp
SimpleArrayOps.hpp
StaticAssert.hpp
string_case_compare.hpp
TeeStreambuf.hpp
TypeList.hpp
  )

APPEND_SET(SOURCES
Bootstrap.cpp
ci_string.cpp
ci_traits.cpp
CSet.cpp
FormatTime.cpp
Identifier.cpp
Marshal.cpp
PrintTable.cpp
  )

IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  stk_util_util
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
