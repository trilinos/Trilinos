#
# Primary Stable options for MPI builds with hybrid GCC 4.5.1 C/C++ and Intel 12.0.4
# Fortran
#

# Included first ot define TRILINOS_TOOLSET_BASE and MKLROOT
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/intel-12.0.8-options.cmake)

SET(INTEL_BASE /opt/intel-12.0.8)

# Point to OpenMPI build with Intel 12.0.8
SET(MPI_BASE_DIR "${INTEL_BASE}/openmpi" CACHE PATH "")
SET(${PROJECT_NAME}_EXTRA_LINK_FLAGS "-Wl,-rpath,${INTEL_BASE}/openmpi/lib ${${PROJECT_NAME}_EXTRA_LINK_FLAGS}" CACHE STRING "")

# Used only when TriKota/Dakota are enabled with CMake
SET(TriKota_ENABLE_DakotaCMake ON CACHE BOOL "")
SET(ENABLE_DAKOTA_TESTS OFF CACHE BOOL "")
SET(HAVE_ACRO OFF CACHE BOOL "")
SET(HAVE_AMPL OFF CACHE BOOL "")
SET(HAVE_X_GRAPHICS OFF CACHE BOOL "")
SET(HAVE_HOPSPACK OFF CACHE BOOL "")
