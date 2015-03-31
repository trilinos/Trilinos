#!/bin/bash

# Simple script I used to drive checkin-test.py from within a CHECKIN subdir of the Trilinos
# source repo, which uses a CHECKIN subdir to store the output files.  This
# makes it easy to do commands like --pull with extra repos to see what the
# status of things are.
#
# Glen Hansen, gahanse@sandia.gov
#
# Here is my directory structure:
#   Codes/
#       Trilinos/
#           packages/
#           checkin-test.py
#           checkin_message
#           build/
#               CHECKIN/
#                   checkin-test-cee-compute011.sh
#               MPI_DEBUG/
#               MPI_RELEASE_SS/
#               SERIAL_RELEASE/
#
# Make sure to "cd ~/Codes/Trilinos/build/CHECKIN; ln -s ../../sampleScripts/checkin-test-cee-compute011.sh ."
#
# To use this script just run it, for example, like (this builds and tests):
#
#  $ ./checkin-test-cee-compute011.sh --do-all
#
# (This pushes after successful tests have been run)
#
#  $ ./checkin-test-cee-compute011.sh --push

for word in "$@"
do
  if [ ${word} == "--help" ]; then
    echo "
To run this script (this builds and tests):

$ ./checkin-test-cee-compute011.sh --do-all

(This pushes after successful tests have been run)

$ ./checkin-test-cee-compute011.sh --push "
exit

  fi
done



#
# Set up configuration files
#
BASE=/projects/albany
BOOSTDIR=$BASE
COMPILER_PATH=/sierra/sntools/SDK/compilers/gcc/4.8.2-RHEL6/bin
MPI_BASE_DIR=/sierra/sntools/SDK/mpi/openmpi/1.6.4-gcc-4.8.2-RHEL6
NETCDF=$BASE
HDFDIR=$BASE
MKL_PATH=/sierra/sntools/SDK/compilers/intel
LABLAS_LIBRARIES="-L$MKL_PATH/lib/intel64 -Wl,--start-group $MKL_PATH/mkl/lib/intel64/libmkl_intel_lp64.a $MKL_PATH/mkl/lib/intel64/libmkl_core.a $MKL_PATH/mkl/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -ldl"

echo "
-D BUILD_SHARED_LIBS:BOOL=OFF
-D Boost_INCLUDE_DIRS:FILEPATH=\"$BOOSTDIR/include\"
-D BoostLib_INCLUDE_DIRS:FILEPATH=\"$BOOSTDIR/include\"
-D Boost_LIBRARY_DIRS:FILEPATH=\"$BOOSTDIR/lib\"
-D BoostLib_LIBRARY_DIRS:FILEPATH=\"$BOOSTDIR/lib\"
-D TPL_ENABLE_Boost:BOOL=ON
-D TPL_ENABLE_BoostLib:BOOL=ON
-D TPL_ENABLE_BLAS:STRING=ON
-D TPL_ENABLE_LAPACK:STRING=ON
-D TPL_BLAS_LIBRARIES:STRING=\"$LABLAS_LIBRARIES\"
-D TPL_LAPACK_LIBRARIES:STRING=\"$LABLAS_LIBRARIES\"
" > COMMON.config

echo "
-DCMAKE_BUILD_TYPE:STRING=DEBUG
-DTPL_ENABLE_MPI:BOOL=ON
-DMPI_BASE_DIR:PATH=$MPI_BASE_DIR
-D TPL_ENABLE_Netcdf:BOOL=ON
-D TPL_ENABLE_HDF5:BOOL=ON
-D Netcdf_INCLUDE_DIRS:PATH=\"$NETCDF/include\"
-D Netcdf_LIBRARY_DIRS:PATH=\"$NETCDF/lib\"
-D HDF5_INCLUDE_DIRS:PATH=\"$HDFDIR/include\"
-D HDF5_LIBRARY_DIRS:PATH=\"$HDFDIR/lib\"
" > MPI_DEBUG.config

#echo "
#-DCMAKE_BUILD_TYPE:STRING=DEBUG
#-DTPL_ENABLE_MPI:BOOL=ON
#-DMPI_BASE_DIR:PATH=$MPI_BASE_DIR
#-D Netcdf_INCLUDE_DIRS:PATH=\"$NETCDF/include\"
#-D Netcdf_LIBRARY_DIRS:PATH=\"$NETCDF/lib\"
#-D HDF5_INCLUDE_DIRS:PATH=\"$HDFDIR/include\"
#-D HDF5_LIBRARY_DIRS:PATH=\"$HDFDIR/lib\"
#" > MPI_DEBUG_SS.config

echo "
-DCMAKE_BUILD_TYPE:STRING=RELEASE
-DTPL_ENABLE_MPI:BOOL=ON
-DMPI_BASE_DIR:PATH=$MPI_BASE_DIR
-DTPL_ENABLE_GLM=OFF
-DTrilinos_ENABLE_Piro:BOOL=OFF
-DTrilinos_ENABLE_Stokhos:BOOL=OFF
-DTrilinos_ENABLE_PyTrilinos:BOOL=OFF
-D Amesos2_ENABLE_KLU2:BOOL=ON
-D TPL_ENABLE_Netcdf:BOOL=ON
-D TPL_ENABLE_HDF5:BOOL=ON
-D TPL_ENABLE_Matio:BOOL=OFF
-D Netcdf_INCLUDE_DIRS:PATH=\"$NETCDF/include\"
-D Netcdf_LIBRARY_DIRS:PATH=\"$NETCDF/lib\"
-D HDF5_INCLUDE_DIRS:PATH=\"$HDFDIR/include\"
-D HDF5_LIBRARY_DIRS:PATH=\"$HDFDIR/lib\"
-D Zlib_INCLUDE_DIRS:PATH=$HDFDIR/include
-D Zlib_LIBRARY_DIRS:PATH=$HDFDIR/lib
-D TPL_HDF5_LIBRARIES:PATH=\"${HDFDIR}/lib/libnetcdf.a;${HDFDIR}/lib/libhdf5_hl.a;${HDFDIR}/lib/libhdf5.a;${HDFDIR}/lib/libz.a\" 
-D Trilinos_EXTRA_LINK_FLAGS:STRING=\"-L${HDFDIR}/lib -lnetcdf -lhdf5_hl -lhdf5 -lz -lgfortran\"
" > MPI_RELEASE_SS.config


echo "
-DCMAKE_BUILD_TYPE:STRING=RELEASE
-DCMAKE_CXX_COMPILER:STRING=$COMPILER_PATH/g++
-DCMAKE_C_COMPILER:STRING=$COMPILER_PATH/gcc
-DCMAKE_Fortran_COMPILER:STRING=$COMPILER_PATH/gfortran
" > SERIAL_RELEASE.config

EXTRA_ARGS=$@

#../../checkin-test.py \
#--make-options="-j 16" \
#--st-extra-builds=MPI_DEBUG_SS \
#--ctest-options="-j 16" \
#--ctest-timeout=1200 \
#--send-email-to="" \
#--no-eg-git-version-check \
#$EXTRA_ARGS

../../checkin-test.py \
--make-options="-j 32" \
--st-extra-builds=MPI_RELEASE_SS \
--ctest-options="-j 16" \
--ctest-timeout=1200 \
--send-email-to="" \
--no-eg-git-version-check \
$EXTRA_ARGS

