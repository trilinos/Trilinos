# Cross-compiler setup
set(CMAKE_SYSTEM_NAME Catamount)
set(CMAKE_SYSTEM_VERSION 1)
set(CMAKE_C_COMPILER cc)
set(CMAKE_C_FLAGS "--target=catamount")
set(CMAKE_CXX_COMPILER CC)
set(CMAKE_CXX_FLAGS "--target=catamount")
set(CMAKE_Fortran_COMPILER ftn)
set(CMAKE_Fortran_FLAGS "--target=catamount")

# Gemini configuration
set(PORTALS_INCLUDE_DIRS /opt/xt-pe/default/include)

# MPI Configuration
set(TPL_ENABLE_MPI ON)
set(MPI_LIBRARY mpich)
set(MPI_EXEC yod)
set(MPI_EXEC_NUMPROCS_FLAG -sz)
set(MPI_USE_COMPILER_WRAPPERS OFF)
set(MPI_INCLUDE_PATH /opt/xt-mpt/default/mpich2-64/T/include)
set(TPL_MPI_INCLUDE_DIRS /opt/xt-mpt/default/mpich2-64/T/include)
