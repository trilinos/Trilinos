# the name of the target operating system
SET(CMAKE_SYSTEM_NAME Catamount)

# set the compiler
#   adding in the -L because apparently the mpi wrapper is missing this path and just adding
#   it to the LD_LIBRARY_PATH doesn't help
set(CMAKE_C_COMPILER mpicc -L/home/sntools/devel/reddish/extras/xt-mpt/2.0.62/sma/P2/lib)
set(CMAKE_CXX_COMPILER mpicxx -L/home/sntools/devel/reddish/extras/xt-mpt/2.0.62/sma/P2/lib)
set(CMAKE_Fortran_COMPILER mpif90 -L/home/sntools/devel/reddish/extras/xt-mpt/2.0.62/sma/P2/lib)

# set the search path for the environment coming with the compiler
# and a directory where you can install your own compiled software
set(CMAKE_FIND_ROOT_PATH
    /home/sntools/devel/reddish/extras/xt-pe/default
    /home/sntools/devel/reddish/extras/xt-mpt/default
    /home/sntools/devel/reddish/extras/xt-mpt/2.0.62/sma/P2
    /home/sntools/devel/reddish/extras/xt-mpt/default/mpich2-64/P2
    /home/sntools/devel/reddish/extras/pgi/7.1.4/linux86-64/7.1-4
  )

# adjust the default behaviour of the FIND_XXX() commands:
# search headers and libraries in the target environment, search
# programs in the host environment
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY BOTH)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE BOTH)
