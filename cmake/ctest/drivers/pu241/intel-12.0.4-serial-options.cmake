#
# Primary Stable options for serial builds with Intel 12.0.4 compiler
#

# Must be including first in order to define TRILINOS_TOOLSET_BASE
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/intel-12.0.4-options.cmake)

SET(TPL_ENABLE_MPI  OFF  CACHE  BOOL "")

SET(CMAKE_CXX_COMPILER     "${INTEL_BIN}/icpc"  CACHE FILEPATH "")
SET(CMAKE_C_COMPILER       "${INTEL_BIN}/icc"   CACHE FILEPATH "")
SET(CMAKE_Fortran_COMPILER "${INTEL_BIN}/ifort" CACHE FILEPATH "")
