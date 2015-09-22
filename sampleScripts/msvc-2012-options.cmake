# Ross Bartlet (2012/12/03): Base options file I use for MS Visual C++ to set
# some basic options that I set with Trilinos_CONFIGURE_OPTIONS_FILE with the
# QT CMake interface.  I use this as a template that I copy into the build
# directory and then I configure with the CMake GUI on Windows by first
# manually setting the option Trilinos_CONFIGURE_OPTIONS_FILE and then running
# configure.  Using this approach makes it easy to reconfigure from scratch.
# This is highly recommended.

SET(CMAKE_BUILD_TYPE DEBUG CACHE STRING "")
SET(BUILD_SHARED_LIBS OFF CACHE BOOL "")
SET(Trilinos_ENABLE_TESTS ON CACHE BOOL "")
SET(Trilinos_ENABLE_Teuchos ON CACHE BOOL "")
SET(Trilinos_ENABLE_RTOp ON CACHE BOOL "")
SET(TPL_ENABLE_Pthread OFF CACHE BOOL "")
