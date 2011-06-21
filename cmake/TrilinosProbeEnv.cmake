INCLUDE(CheckIncludeFileCXX)

IF(MSVC)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE 
    -D_CRT_NONSTDC_NO_DEPRECATE  -D_SCL_SECURE_NO_WARNINGS)
  INCLUDE_DIRECTORIES(${Trilinos_SOURCE_DIR}/commonTools/WinInterface/include)
  # find the CLAPACK built by CMake on the machine for MSVC
  # if found it will set the BLAS and LAPACK libraries
  FIND_PACKAGE(CLAPACK 3.2.1 NO_MODULE)
  IF(CLAPACK_FOUND)
    SET(TPL_BLAS_LIBRARIES blas CACHE INTERNAL "")
    SET(TPL_LAPACK_LIBRARIES lapack CACHE INTERNAL "")
  ENDIF()
ENDIF()

IF (WIN32 AND NOT CYGWIN)
  SET(NATIVE_MS_WINDOWS TRUE)
ELSE()
  SET(NATIVE_MS_WINDOWS FALSE)
ENDIF()

# Probe for non-standard headers

IF (Trilinos_ENABLE_CXX)
  CHECK_INCLUDE_FILE_CXX(sys/time.h HAVE_SYS_TIME_H)
  CHECK_INCLUDE_FILE_CXX(time.h HAVE_TIME_H)
  CHECK_INCLUDE_FILE_CXX(stdint.h HAVE_STDINT_H)
  CHECK_INCLUDE_FILE_CXX(inttypes.h HAVE_INTTYPES_H)
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

# Determine C++-0x supported features

IF (Trilinos_ENABLE_CXX11)
  INCLUDE(TrilinosCXX11Support)
  CHECK_CXX11_SUPPORT(Trilinos_ENABLE_CXX11)
  MESSAGE("-- Trilinos_ENABLE_CXX11=${Trilinos_ENABLE_CXX11}")
ENDIF()

# Set up some MPI info

IF (TPL_ENABLE_MPI)
  SET(HAVE_MPI TRUE)
ELSE()
  SET(HAVE_MPI FALSE)
ENDIF()

# OpenMP isn't really a TPL because support is built into the compiler.

IF(Trilinos_ENABLE_OpenMP)
  INCLUDE(FindOpenMP)
  IF(OPENMP_FOUND)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    # FindOpenMP.cmake doesn't find Fortran flags.  Mike H said this is safe.
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
  ELSE()
    MESSAGE(FATAL_ERROR "Could not find OpenMP, try setting OpenMP_C_FLAGS and OpenMP_CXX_FLAGS directly")
  ENDIF(OPENMP_FOUND)
ENDIF(Trilinos_ENABLE_OpenMP)

# Check if we need the math library or not and find the right one
IF (NOT NATIVE_MS_WINDOWS)
  INCLUDE(MathLibraryNeeded)
ENDIF()

# Check for isnan and isinf support
IF (${PROJECT_NAME}_ENABLE_CXX)
  INCLUDE(FiniteValue)
ENDIF()

# Check for Doxygen/dot - We can use variables set in this check to
# enable/disable the grapical dependency graphs in doxygen Doxyfiles.
INCLUDE(FindDoxygen)

