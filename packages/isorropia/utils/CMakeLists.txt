
INCLUDE(TribitsLibraryMacros)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS)
SET(SOURCES)

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS ${HEADERS}
  ispatest_epetra_utils.hpp
  ispatest_lbeval_utils.hpp
  ispatest_utils.hpp
  ispatest_read_distribution.hpp
  )

SET(SOURCES ${SOURCES}
  ispatest_epetra_utils.cpp
  ispatest_lbeval_utils.cpp
  ispatest_utils.cpp
  ispatest_read_distribution.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  isorropia_test_utils
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES}
  TESTONLY
  DEPLIBS isorropia
  )
