
INCLUDE(TrilinosPackageLibraryHelpers)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS)
SET(SOURCES)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

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

TRILINOS_PACKAGE_ADD_LIBRARY(
  isorropia_test_utils
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS isorropia
  )
# 2008/11/20: rabartl: ToDo: Above: Add an option like TESTLIB to the
# above macro call to in order to result in the library not being added
# to the main libraries and so that it will not be installed!

#
# D) Export the dependancy variables of this package for other packages to use
#

TRILINOS_PACKAGE_EXPORT_DEPENDENCY_VARIABLES()
