# This file contains the options needed to both run the pull request testing
# for Trilinos for the Linux Intel 19 pull request testing builds, and to reproduce
# the errors reported by those builds. Prior to using this this file, the
# appropriate set of SEMS modules must be loaded and accessible through the
# SEMS NFS mount. (See the sems/PullRequestGCC*TestingEnv.sh files.)

# Usage: cmake -C PullRequestLinuxIntel19.0.5TestingSettings.cmake

set (CMAKE_CXX_STANDARD "14" CACHE STRING "Set C++ standard to C++14")

include("${CMAKE_CURRENT_LIST_DIR}/PullRequestLinuxCommonTestingSettings.cmake")

set (Tpetra_INST_INT_INT ON CACHE BOOL "INST_INT_INT ON")
set (Trilinos_ENABLE_STKBalance OFF CACHE BOOL "Hard disabled since Tpetra_INST_INT_INT=ON in this build" FORCE)
#STK-TODO: try to remember to come back and remove this when stk-balance
#is able to tolerate int as a global-index.

set (TPL_Netcdf_LIBRARIES "-L${Netcdf_LIBRARY_DIRS}/lib;${Netcdf_LIBRARY_DIRS}/libnetcdf.so;${Netcdf_LIBRARY_DIRS}/libpnetcdf.a" CACHE STRING "Set by default for CUDA PR testing")

set(CMAKE_CXX_FLAGS "-fPIC -Wall -Warray-bounds -Wchar-subscripts -Wcomment -Wenum-compare -Wformat -Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull -Wparentheses -Wpointer-sign -Wreorder -Wreturn-type -Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-function -Wunused-but-set-variable -Wunused-variable -Wwrite-strings" CACHE STRING "enable relocatable code, turn on WError settings")

# -lifcore is needed for CMake > 3.10.x builds.
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lifcore" CACHE STRING "updated by Pull Request")

#set (Anasazi_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Belos_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Domi_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Epetra_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (EpetraExt_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (FEI_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Ifpack_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Ifpack2_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Intrepid2_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Kokkos_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (KokkosKernels_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (ML_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (MueLu_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (NOX_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Panzer_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Phalanx_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Pike_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Piro_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (ROL_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Sacado_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (SEACAS_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Shards_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Stokhos_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Tempus_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Teuchos_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Tpetra_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Triutils_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Zoltan2_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
