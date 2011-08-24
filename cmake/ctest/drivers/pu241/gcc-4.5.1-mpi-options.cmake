# Primary Stable options for MPI builds with GCC 4.5.1
INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/gcc-4.5.1-base-options.cmake)
# this path refers to the pure GNU OpenMPI installation
SET(MPI_BASE_DIR "${TRILINOS_TOOLSET_BASE}" CACHE PATH "" FORCE)
# this path refers to the GNU C/C++ and Intel Fortran OpenMPI installation
# SET(MPI_BASE_DIR "${TRILINOS_TOOLSET_BASE}/ifort-12.0.4" CACHE PATH "" FORCE)
# SET(Trilinos_EXTRA_LINK_FLAGS "-Wl,-rpath,${TRILINOS_TOOLSET_BASE}/ifort-12.0.4/lib64 ${Trilinos_EXTRA_LINK_FLAGS}" CACHE STRING "")

# Used only when TriKota/Dakota are enabled with CMake
SET(TriKota_ENABLE_DakotaCMake ON CACHE BOOL "" FORCE)
SET(ENABLE_DAKOTA_TESTS OFF CACHE BOOL "" FORCE)
SET(HAVE_ACRO OFF CACHE BOOL "" FORCE)
SET(HAVE_AMPL OFF CACHE BOOL "" FORCE)
SET(HAVE_X_GRAPHICS OFF CACHE BOOL "" FORCE)
SET(HAVE_HOPSPACK OFF CACHE BOOL "" FORCE)
