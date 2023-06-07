rm -fr CMake*

module purge

INSTALL_DIR=/ascldap/users/knliege/local/trilinos_all/pytrilinos_test
TRILINOS_DIR=/ascldap/users/knliege/dev/trilinos_all/Trilinos

export ATDM_CONFIG_REGISTER_CUSTOM_CONFIG_DIR=${TRILINOS_DIR}/cmake/std/atdm/contributed/weaver
source $TRILINOS_DIR/cmake/std/atdm/load-env.sh weaver-cuda-10.1-opt

cmake \
-D LAPACK_LIBRARY_DIRS=${LAPACK_ROOT} \
-D BLAS_LIBRARY_DIRS=${LAPACK_ROOT} \
-D Trilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake \
-D CMAKE_Fortran_COMPILER:FILEPATH=mpif77 \
-D BUILD_SHARED_LIBS:BOOL=${ATDM_CONFIG_SHARED_LIBS} \
-D CMAKE_INSTALL_PREFIX:FILEPATH=${INSTALL_DIR} \
-D CMAKE_BUILD_TYPE:STRING=RELEASE \
-D TPL_ENABLE_CUDA:BOOL=ON \
-D Trilinos_ENABLE_PyTrilinos2:BOOL=ON \
-D PyTrilinos2_BINDER_EXECUTABLE=/ascldap/users/knliege/dev/binder/prefix/build/bin/binder \
-D PyTrilinos2_BINDER_GCC_TOOLCHAIN=/home/projects/ppc64le/gcc/7.2.0 \
-D PyTrilinos2_BINDER_FLAGS="--cuda-gpu-arch=sm_70;--cuda-path=/home/projects/ppc64le-pwr9-nvidia/cuda/10.1.105" \
-D PyTrilinos2_ENABLE_TESTS=ON \
-D PyTrilinos2_ENABLE_Binder=ON \
-D PyTrilinos2_ENABLE_Update_Binder=OFF \
-D PYTHON_EXECUTABLE="~/weaver/miniconda3/bin/python" \
-D Trilinos_ENABLE_Tpetra:BOOL=ON \
-D Trilinos_ENABLE_Epetra:BOOL=OFF \
-D BUILD_SHARED_LIBS:BOOL=ON \
-D Tpetra_ENABLE_TESTS=ON \
-D Amesos2_ENABLE_TESTS=ON \
-D Trilinos_ENABLE_Fortran:BOOL=OFF \
-D Trilinos_EXTRA_LINK_FLAGS:STRING="-lgfortran" \
-D Kokkos_ENABLE_CUDA_UVM=OFF \
-D Trilinos_ENABLE_OpenMP=OFF \
-D Kokkos_ENABLE_OPENMP:BOOL=OFF \
\
$TRILINOS_DIR
