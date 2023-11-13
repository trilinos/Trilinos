################################################################################
#
# Set up cts1 env for EMPIRE ATDM builds of Trilinos
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
. /projects/sems/modulefiles/utils/sems-archive-modules-init.sh
module load sems-archive-env
module load sems-archive-ninja_fortran/1.8.2
module load sems-archive-git/2.10.1
module load sems-archive-cmake/3.19.1

if [ "$ATDM_CONFIG_NODE_TYPE" == "OPENMP" ] ; then
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=8
  export OMP_NUM_THREADS=2
else
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=16
fi

sparc_tpl_base=${ATDM_CONFIG_SPARC_TPL_BASE}

if [ "$ATDM_CONFIG_COMPILER" == "INTEL-18.0.2_OPENMPI-4.0.1" ]; then
    module load intel/18.0.2.199
    module load mkl/18.0.5.274 # Needed to address defect in MKL (#5316, #3992, #3914)
    module load openmpi-intel/4.0

    sparc_tpl_ext=cts1-bdw_intel-19.0.5
    sparc_tpl_mpi_ext=cts1-bdw_intel-19.0.5_openmpi-4.0.1

    export PATH=/usr/tce/packages/gcc/gcc-6.1.0/bin:${PATH}
    export LD_LIBRARY_PATH=/usr/tce/packages/gcc/gcc-6.1.0/lib64:${LD_LIBRARY_PATH}
else
    echo
    echo "***"
    echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not supported on this system!"
    echo "***"
    return
fi

export OMPI_CXX=`which icpc`
export OMPI_CC=`which icc`
export OMPI_FC=`which ifort`
export ATDM_CONFIG_LAPACK_LIBS="-mkl"
export ATDM_CONFIG_BLAS_LIBS="-mkl"

export BOOST_ROOT=${sparc_tpl_base}/boost-1.72.0/00000000/${sparc_tpl_ext}
export HDF5_ROOT=${sparc_tpl_base}/hdf5-1.10.5/00000000/${sparc_tpl_mpi_ext}
export CGNS_ROOT=${sparc_tpl_base}/cgns-c09a5cd/d313cc2f822078e47c7dbdee074ecb0431e573eb/${sparc_tpl_mpi_ext}
export PNETCDF_ROOT=${sparc_tpl_base}/pnetcdf-1.12.1/6144dc67b2041e4093063a04e89fc1e33398bd09/${sparc_tpl_mpi_ext}
export NETCDF_ROOT=${sparc_tpl_base}/netcdf-4.7.0/24baa07a3fa1ff9dbc8e70dc591ebbdec56783b2/${sparc_tpl_mpi_ext}
export PARMETIS_ROOT=${sparc_tpl_base}/parmetis-4.0.3/00000000/${sparc_tpl_mpi_ext}
export METIS_ROOT=${sparc_tpl_base}/parmetis-4.0.3/00000000/${sparc_tpl_mpi_ext}
export LIBHIO_ROOT=${sparc_tpl_base}/libhio-1.4.1.2/00000000/${sparc_tpl_mpi_ext}
export EUCLID_ROOT=${sparc_tpl_base}/euclid-20.23/8b68b12f72b59648c9a0a962a6d55ea978199860/${sparc_tpl_mpi_ext}
export SGM_ROOT=${sparc_tpl_base}/sgm-20.23/00000000/${sparc_tpl_mpi_ext}
export SUPERLUDIST_ROOT=${sparc_tpl_base}/superlu_dist-5.4.0/a3121eaff44f7bf7d44e625c3b3d2a9911e58876/${sparc_tpl_mpi_ext}

export ATDM_CONFIG_USE_HWLOC=OFF
export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"
export ATDM_CONFIG_NETCDF_LIBS="-L${NETCDF_ROOT}/lib64;-L${PNETCDF_ROOT}/lib;${NETCDF_ROOT}/lib64/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${ATDM_CONFIG_HDF5_LIBS};-lcurl"
export ATDM_CONFIG_SUPERLUDIST_INCLUDE_DIRS=${SUPERLUDIST_ROOT}/include
export ATDM_CONFIG_SUPERLUDIST_LIBS=${SUPERLUDIST_ROOT}/lib64/libsuperlu_dist.a
export ATDM_CONFIG_BINUTILS_LIBS="/usr/lib64/libbfd.so;/usr/lib64/libiberty.a"

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
export FC=mpif77
export F90=mpif90

# Define function atdm_run_script_on_compute_node
source $ATDM_SCRIPT_DIR/common/define_run_on_slurm_compute_node_func.sh

# Point CMake 3.18 compiler checks to missing symbols
export LDFLAGS="$LDFLAGS -lifcore"

export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
