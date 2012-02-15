# Cross-compiler setup
set(CMAKE_C_COMPILER cc)
set(CMAKE_CXX_COMPILER CC)
set(CMAKE_Fortran_COMPILER ftn)

# Gemini configuration
#set(GEMINI_INCLUDE_DIRS /opt/cray/gni-headers/default/include)
set(GEMINI_INCLUDE_DIRS /opt/cray/gni-headers/2.1-1.0400.4156.6.1.gem/include)

# MPI Configuration
set(TPL_ENABLE_MPI ON)
set(MPI_LIBRARY mpich)
set(MPI_EXEC aprun)
set(MPI_EXEC_NUMPROCS_FLAG -n)
set(MPI_USE_COMPILER_WRAPPERS OFF)
#set(MPI_INCLUDE_PATH /opt/cray/mpt/default/xt/gemini/mpich2-pgi/include)
#set(TPL_MPI_INCLUDE_DIRS /opt/cray/mpt/default/xt/gemini/mpich2-pgi/include)
set(MPI_INCLUDE_PATH /opt/cray/mpt/5.4.2/xt/gemini/mpich2-pgi/include)
set(TPL_MPI_INCLUDE_DIRS /opt/cray/mpt/5.4.2/xt/gemini/mpich2-pgi/include)
