#
# Base options for all GCC 4.8.3 builds
#

# Define the core compilers
IF (NOT TRILINOS_TOOLSET_BASE)
  SET(TRILINOS_TOOLSET_BASE  /projects/vera/gcc-4.8.3/toolset)
ENDIF()
SET(GCC_BASE_DIR ${TRILINOS_TOOLSET_BASE}/gcc-4.8.3)

# Point to the right MPI
SET(MPI_BASE_DIR "${TRILINOS_TOOLSET_BASE}/mpich-3.1.3" CACHE PATH "")

IF (NOT TPL_ENABLE_MPI)
  SET(CMAKE_C_COMPILER gcc)
  SET(CMAKE_CXX_COMPILER g++)
  SET(CMAKE_Fortran_COMPILER gfortran)
ENDIF()

# Add rpath for compiler libraries and gomp for parts built with OpenMP
SET(${PROJECT_NAME}_EXTRA_LINK_FLAGS
  "-lgomp -Wl,-rpath,${GCC_BASE_DIR}/lib64"
  CACHE STRING "")

# Build shared libs by default to save massive amounts of disk space
SET(BUILD_SHARED_LIBS  ON  CACHE  BOOL
  "Set by default in gcc-4.8.3-base-options.cmake")

# Turn on explicit template instantaition by default
SET(Trilinos_ENABLE_EXPLICIT_INSTANTIATION  ON  CACHE BOOL
  "Set by default in gcc-4.8.3-base-options.cmake")

# Turn on KLU2 in Amesos2
SET(Amesos2_ENABLE_KLU2  ON  CACHE BOOL
  "Set by default in gcc-4.8.3-base-options.cmake")

# Set up valgrind options
SET( MEMORYCHECK_COMMAND
  /projects/vera/common_tools/valgrind-3.9.0/bin/valgrind
  CACHE  FILEPATH
  "Set by default in gcc-4.8.3-base-options.cmake")
SET( MEMORYCHECK_COMMAND_OPTIONS
  "-q --trace-children=yes --tool=memcheck --leak-check=yes --leak-check=full --workaround-gcc296-bugs=yes --num-callers=50"
  CACHE  STRING
  "Set by default in gcc-4.8.3-base-options.cmake")
