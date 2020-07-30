TRIBITS_PACKAGE(Gtest)

# Option to skip installation of gtest package
SET(${PACKAGE_NAME}_SKIP_INSTALL FALSE CACHE BOOL
  "Set to TRUE to skip the install of Gtest header and library")
IF (${PACKAGE_NAME}_SKIP_INSTALL)
  SET(Gtest_NO_INSTALL_LIB_OR_HEADERS_ARG NO_INSTALL_LIB_OR_HEADERS CACHE INTERNAL "")
  SET(INSTALL_GTEST CACHE BOOL OFF) # referenced by internal googletest package
ELSE()
  SET(Gtest_NO_INSTALL_LIB_OR_HEADERS_ARG CACHE INTERNAL "")
ENDIF()

# Configure googletest
add_subdirectory(googletest)

# Create dummy tribits library linking to internal gtest library so tribits
# will set dependencies of downstream packages correctly.  Unfortunately we
# appear to have tp provide some source file
TRIBITS_ADD_LIBRARY(
  trilinos-gtest
  SOURCES trilinos-gtest.cpp
  ${Gtest_NO_INSTALL_LIB_OR_HEADERS_ARG}
  ADDED_LIB_TARGET_NAME_OUT TrilinosGtestLibName
  )
TARGET_LINK_LIBRARIES("${TrilinosGtestLibName}" gtest)

TRIBITS_PACKAGE_POSTPROCESS()
