INCLUDE(TribitsLibraryMacros)
#INCLUDE(TribitsAddExecutableAndTest)

SET(HEADERS "")
SET(SOURCES "")

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  Tempus_ConvergenceTestUtils.hpp
  )

APPEND_SET(SOURCES
  )

# This is a header only library at the moment. Some versions of ar
# fail to create the library with an error if it is empty. So we are
# commenting out this library declaration until we have a true object
# defined in it. If we always used ETI, then we could uncomment this.

#TRIBITS_ADD_LIBRARY(
#  tempus_test_utils
#  HEADERS ${HEADERS}
#  SOURCES ${SOURCES}
#  TESTONLY
#  )
