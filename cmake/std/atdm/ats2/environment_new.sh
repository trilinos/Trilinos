################################################################################
#
# Set up env on ats2 (vortex) for ATMD builds of Trilinos
#
# This source script gets the settings from the ATDM_CONFIG_BUILD_NAME var.
#
################################################################################

echo "Using $ATDM_CONFIG_SYSTEM_NAME toss3 compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

export ATDM_CONFIG_ENABLE_SPARC_SETTINGS=ON
export ATDM_CONFIG_USE_NINJA=ON
export ATDM_CONFIG_BUILD_COUNT=8
# export ATDM_CONFIG_CMAKE_JOB_POOL_LINK=2
# NOTE: Above, currently setting CMAKE_JOB_POOL_LINK results in a build
# failures with Ninja.  See https://gitlab.kitware.com/snl/project-1/issues/60

# We do this twice since sems modules are wacked and we get errors to the screen on a purge
# The second purge will catch any real errors with purging ...
module purge &> /dev/null
module purge
#. /projects/sems/modulefiles/utils/sems-modules-init.sh
#module load sems-env
#module load atdm-env
module load git/2.20.0

if [ "$ATDM_CONFIG_NODE_TYPE" == "OPENMP" ] ; then
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=8
  export OMP_NUM_THREADS=2
else
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=16
fi

# Common paths and modules for both intel-1{8,9}
sparc_tpl_base=${ATDM_CONFIG_SPARC_TPL_BASE}

if [ "$ATDM_CONFIG_COMPILER" == "GNU-7.3.1_SPMPI-2019.06.24" ]; then
  module load cmake/3.14.5 # Vortex does not have 3.12.2
  module load gcc/7.3.1 # Make spectrum-mpi available
#    module load mkl/18.0.5.274 # Needed to find libmkl_intel_lp64.so (ATDV-212)
    #:module show sparc-dev/gcc-7.3.1_spmpi-2019.06.24:
    #module load sparc-tools/aerotools; module unload xl; module unload spectrum-mpi
    #module load cmake/3.14.5
    #export PATH=/projects/sparc/tools/ats2/clang-8.0.1/bin:$PATH
    module load spectrum-mpi/2019.06.24
    #...
    module load lapack/3.8.0-gcc-4.9.3

    sparc_tpl_ext=ats2-pwr9_gcc-7.3.1
    sparc_tpl_mpi_ext=ats2-pwr9_gcc-7.3.1_spmpi-2019.06.24
    # rabartl: ToDo: Above, we need to find a way to extract 'cts1-bdw' out of
    # this file for this to be general!

    export CBLAS_ROOT=/usr/tcetmp/packages/lapack/lapack-3.8.0-gcc-4.9.3 #
    export LAPACK_ROOT=/usr/tcetmp/packages/lapack/lapack-3.8.0-gcc-4.9.3 #
    export COMPILER_ROOT=/usr/tce/packages/gcc/gcc-7.3.1
    export SPARC_HDF5=hdf5-1.10.5

    export PATH=${COMPILER_ROOT}/bin:${PATH}
    export LD_LIBRARY_PATH=${COMPILER_ROOT}/lib:${LD_LIBRARY_PATH}
    export BINUTILS_ROOT=${COMPILER_ROOT}
    export LIBRARY_PATH=${BINUTILS_ROOT}/lib:${LIBRARY_PATH}
    export LIBRARY_PATH=${CBLAS_ROOT}/lib:${LIBRARY_PATH}
    export INCLUDE=${BINUTILS_ROOT}/include:${INCLUDE}
    export CPATH=${BINUTILS_ROOT}/include:${CPATH}
    export ATDM_CONFIG_BINUTILS_LIBS="${BINUTILS_ROOT}/lib/libbfd.a;${BINUTILS_ROOT}/lib/libiberty.a"

    #export CC=mpicc
    #export CXX=mpicxx
    export F77=mpifort
    export FC=mpifort
    export F90=mpifort
else
    echo
    echo "***"
    echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not supported on this system!"
    echo "***"
    return
fi

#TODO: don't use intel for this... do we need a new environment_new.sh file?
#export ATDM_CONFIG_LAPACK_LIBS="-mkl"
#export ATDM_CONFIG_BLAS_LIBS="-mkl"
export ATDM_CONFIG_LAPACK_LIBS="-L${LAPACK_ROOT}/lib;-llapack;-lgfortran;-lgomp"
export ATDM_CONFIG_BLAS_LIBS="-L${BLAS_ROOT}/lib;-lblas;-lgfortran;-lgomp;-lm"

export BOOST_ROOT=${sparc_tpl_base}/boost-1.65.1/00000000/${sparc_tpl_ext}
export HDF5_ROOT=${sparc_tpl_base}/hdf5-1.10.5/00000000/${sparc_tpl_mpi_ext}
export CGNS_ROOT=${sparc_tpl_base}/cgns-c09a5cd/27e5681f1b74c679b5dcb337ac71036d16c47977/${sparc_tpl_mpi_ext}
export PNETCDF_ROOT=${sparc_tpl_base}/pnetcdf-1.10.0/6144dc67b2041e4093063a04e89fc1e33398bd09/${sparc_tpl_mpi_ext}
export NETCDF_ROOT=${sparc_tpl_base}/netcdf-4.7.0/58bc48d95be2cc9272a18488fea52e1be1f0b42a/${sparc_tpl_mpi_ext}
export PARMETIS_ROOT=${sparc_tpl_base}/parmetis-4.0.3/00000000/${sparc_tpl_mpi_ext}
export METIS_ROOT=${sparc_tpl_base}/parmetis-4.0.3/00000000/${sparc_tpl_mpi_ext}
#export LIBHIO_ROOT=${sparc_tpl_base}/libhio-1.4.1.2/00000000/${sparc_tpl_mpi_ext}
#export EUCLID_ROOT=${sparc_tpl_base}/euclid-19.30/95b8242ad729449fdb86591b9643c11463b61fa2/${sparc_tpl_mpi_ext}
#export SGM_ROOT=${sparc_tpl_base}/sgm-19.30/00000000/${sparc_tpl_mpi_ext}
export SUPERLUDIST_ROOT=${sparc_tpl_base}/superlu_dist-5.4.0/a3121eaff44f7bf7d44e625c3b3d2a9911e58876/${sparc_tpl_mpi_ext}

export ATDM_CONFIG_USE_HWLOC=OFF
export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"
export ATDM_CONFIG_NETCDF_LIBS="-L${BOOST_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${PNETCDF_ROOT}/lib;-L${HDF5_ROOT}/lib;${BOOST_ROOT}/lib/libboost_program_options.a;${BOOST_ROOT}/lib/libboost_system.a;${NETCDF_ROOT}/lib/libnetcdf.a;${NETCDF_ROOT}/lib/libpnetcdf.a;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl;-lcurl"
export ATDM_CONFIG_SUPERLUDIST_INCLUDE_DIRS=${SUPERLUDIST_ROOT}/include
export ATDM_CONFIG_SUPERLUDIST_LIBS=${SUPERLUDIST_ROOT}/lib64/libsuperlu_dist.a

# not sure what below does.  It was in the original environment script
#unset ATTB_ENV

# Set MPI wrappers
export MPICC=`which mpicc`
export MPICXX=`which mpicxx`
export MPIF90=`which mpif90`

export ATDM_CONFIG_MPI_EXEC=srun
export ATDM_CONFIG_MPI_PRE_FLAGS="--mpi=pmi2;--ntasks-per-node;36"
export ATDM_CONFIG_MPI_EXEC_NUMPROCS_FLAG=--ntasks

# Set the default compilers
export CC=mpicc
export CXX=mpicxx

# Define function atdm_run_script_on_compute_node
source $ATDM_SCRIPT_DIR/common/define_run_on_slurm_compute_node_func.sh

export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
