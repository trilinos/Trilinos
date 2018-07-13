August 10, 2017
---------------

At present, I am uncertain as how to automatically find Eigen using the usual 
TRIBITS_TPL_INCLUDE_DIRS_AND_LIBRARIES call. My suspicion is that this is because the 
typical include directives are for Eigen/Core, Eigen/Dense, etc, which are files that 
contain relative include paths for all of the headers in the directories
../Eigen/src/Core and ../Eigen/src/Dense respectively. On the other hand, the 
cmake function find_package(Eigen3 REQUIRED) succeeds in setting EIGEN3_INCLUDE_DIR
(typically /usr/lib/eigen3 on Debian flavored systems), although this does not conform
to the standard Trilinos approach to TPLs. 

For the time being, I have found that configuring Trilinos with a script like the
following seems to work:

    CXXFLAGS="-std=c++11 -fPIC"
  
    cmake \
     -D TPL_ENABLE_LIBRARIES_AND_HEADERS:BOOL=ON              \
       -D TPL_BLAS_LIBRARIES:FILEPATH='/usr/lib/libblas.so.3' \
       -D TPL_LAPACK_LIBRARIES:PATH='/usr/lib/liblapack.so.3' \
     -D TPL_ENABLE_MPI:BOOL=ON                                \
     -D TPL_ENABLE_Eigen:BOOL=ON                              \
       -D Eigen_INCLUDE_DIRS:PATH='/usr/include/eigen3'       \
     -D CMAKE_BUILD_TYPE:STRING=DEBUG                         \
     -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF        \
     -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=OFF       \
     -D Trilinos_ENABLE_TESTS:BOOL=OFF                        \
     -D Trilinos_ENABLE_ROL:BOOL=ON                           \
       -D ROL_ENABLE_TESTS:BOOL=ON                            \
    ${TRILINOS_HOME}


