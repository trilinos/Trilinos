INCLUDE(CheckIncludeFileCXX)

INCLUDE(CMakeBuildTypesList)

IF (WIN32 AND NOT CYGWIN)
  SET(NATIVE_MS_WINDOWS TRUE)
ELSE()
  SET(NATIVE_MS_WINDOWS FALSE)
ENDIF()

# Set to release build by default

IF (NOT CMAKE_BUILD_TYPE)
  MESSAGE(STATUS "Setting CMAKE_BUILD_TYPE=RELEASE since none was set")
  SET(CMAKE_BUILD_TYPE RELEASE)
ELSE()
  LIST(FIND CMAKE_BUILD_TYPES_LIST ${CMAKE_BUILD_TYPE} BUILD_TYPE_IDX)
  IF (BUILD_TYPE_IDX EQUAL -1)
    MESSAGE(SEND_ERROR "Error, the given CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
      " is not in the list of valid values \"${CMAKE_BUILD_TYPES_LIST}\"!")
  ENDIF()
ENDIF()
PRINT_VAR(CMAKE_BUILD_TYPE)

# Probe for non-standard headers

IF (Trilinos_ENABLE_CXX)
  CHECK_INCLUDE_FILE_CXX(sys/time.h HAVE_SYS_TIME_H)
  CHECK_INCLUDE_FILE_CXX(time.h HAVE_TIME_H)
ENDIF()

SET(HAVE_ALGORITHM TRUE)
SET(HAVE_CASSERT TRUE)
SET(HAVE_CCTYPE TRUE)
SET(HAVE_CERRNO TRUE)
SET(HAVE_CLIMITS TRUE)
SET(HAVE_CMATH TRUE)
SET(HAVE_COMPLEX TRUE)
SET(HAVE_CSTDARG TRUE)
SET(HAVE_CSTDIO TRUE)
SET(HAVE_CSTDLIB TRUE)
SET(HAVE_CSTRING TRUE)
SET(HAVE_IOMANIP TRUE)
SET(HAVE_IOSTREAM TRUE)
SET(HAVE_ITERATOR TRUE)
SET(HAVE_LIST TRUE)
SET(HAVE_MAP TRUE)
SET(HAVE_MEMORY TRUE)
SET(HAVE_MUTABLE TRUE)
SET(HAVE_NAMESPACES TRUE)
SET(HAVE_NEW_FOR_SCOPING TRUE)
SET(HAVE_NUMERIC TRUE)
SET(HAVE_NUMERIC_LIMITS TRUE)
SET(HAVE_POW TRUE)
SET(HAVE_SET TRUE)
SET(HAVE_SSTREAM TRUE)
SET(HAVE_FSTREAM TRUE)
SET(HAVE_STDEXCEPT TRUE)
SET(HAVE_STRING TRUE)
SET(HAVE_VECTOR TRUE)

# 2008/12/20: rabartl: Above: All of these defines should be removed
# because we decided that we were going to assume that all compilers
# have these C++98 standard features.  We will deal with cases where
# this is not true but we should not assume the worst right from the
# beginning.

# Find Perl

FIND_PACKAGE(Perl)

# Do Fortran stuff

INCLUDE(TrilinosFortranMangling)

# Get BLAS name mangling
 
INCLUDE(TrilinosBLASMangling)

# Set up some MPI info

IF (TPL_ENABLE_MPI)
  SET(HAVE_MPI TRUE)
  ADD_DEFINITIONS(-DMPICH_IGNORE_CXX_SEEK)
ELSE()
  SET(HAVE_MPI FALSE)
ENDIF()

# Check if we need the math library or not and find the right one
IF (NOT NATIVE_MS_WINDOWS)
  INCLUDE(MathLibraryNeeded)
ENDIF()

# Check for isnan and isinf support
INCLUDE(FiniteValue)

# Check for Doxygen/dot - We can use variables set in this check to
# enable/disable the grapical dependency graphs in doxygen Doxyfiles.
INCLUDE(FindDoxygen)
