################################################################################
#
# Set up env on toss3 (chama and serrano) for ATMD builds of Trilinos
#
# This source script gets the settings from the ATDM_CONFIG_BUILD_NAME var.
#
################################################################################

if [ "$ATDM_CONFIG_COMPILER" == "DEFAULT" ] ; then
  export ATDM_CONFIG_COMPILER=INTEL
fi

echo "Using toss3 compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

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
. /projects/sems/modulefiles/utils/sems-modules-init.sh
module load sems-env
module load atdm-env
module load atdm-ninja_fortran/1.7.2

if [ "$ATDM_CONFIG_NODE_TYPE" == "OPENMP" ] ; then
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=8
  export OMP_NUM_THREADS=2
else
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=16
fi

# Common paths and modules for both intel-1{8,9}
sparc_tpl_base=/projects/sparc/tpls/cts1-bdw
module load mkl/18.0.5.274
module load cmake/3.12.2

if [ "$ATDM_CONFIG_COMPILER" == "INTEL-18.0.2_OPENMPI-2.0.3" ]; then
    module load intel/18.0.2.199
    module load openmpi-intel/2.0

    sparc_tpl_ext=cts1-bdw_intel-18.0.2
    sparc_tpl_mpi_ext=cts1-bdw_intel-18.0.2_openmpi-2.0.3

    export PATH=/usr/tce/packages/gcc/gcc-6.1.0/bin:${PATH}
    export LD_LIBRARY_PATH=/usr/tce/packages/gcc/gcc-6.1.0/lib64:${LD_LIBRARY_PATH}
elif [ "$ATDM_CONFIG_COMPILER" == "INTEL-19.0.5_OPENMPI-4.0.1" ]; then
    module load intel/19.0.5.281
    module load openmpi-intel/4.0

    sparc_tpl_ext=cts1-bdw_intel-19.0.5
    sparc_tpl_mpi_ext=cts1-bdw_intel-19.0.5_openmpi-4.0.1

    export CBLAS_ROOT=/projects/global/toss3/compilers/intel/intel_2019/compilers_and_libraries_2019.5.281/linux
    export COMPILER_ROOT=/projects/global/toss3/compilers/intel/intel_2017/compilers_and_libraries_2019.5.281/linux
    export SPARC_HDF5=hdf5-1.10.5

    export PATH=/usr/tce/packages/gcc/gcc-4.9.3/bin:${PATH}
    export LD_LIBRARY_PATH=/usr/tce/packages/gcc/gcc-4.9.3/lib64:${LD_LIBRARY_PATH}

    export F77=mpif77
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

export BOOST_ROOT=${sparc_tpl_base}/boost-1.65.1/00000000/${sparc_tpl_ext}
export HDF5_ROOT=${sparc_tpl_base}/hdf5-1.10.5/00000000/${sparc_tpl_mpi_ext}
export CGNS_ROOT=${sparc_tpl_base}/cgns-c09a5cd/27e5681f1b74c679b5dcb337ac71036d16c47977/${sparc_tpl_mpi_ext}
export PNETCDF_ROOT=${sparc_tpl_base}/pnetcdf-1.10.0/6144dc67b2041e4093063a04e89fc1e33398bd09/${sparc_tpl_mpi_ext}
export NETCDF_ROOT=${sparc_tpl_base}/netcdf-4.7.0/58bc48d95be2cc9272a18488fea52e1be1f0b42a/${sparc_tpl_mpi_ext}
export PARMETIS_ROOT=${sparc_tpl_base}/parmetis-4.0.3/00000000/${sparc_tpl_mpi_ext}
export METIS_ROOT=${sparc_tpl_base}/parmetis-4.0.3/00000000/${sparc_tpl_mpi_ext}
export LIBHIO_ROOT=${sparc_tpl_base}/libhio-1.4.1.2/00000000/${sparc_tpl_mpi_ext}
export EUCLID_ROOT=${sparc_tpl_base}/euclid-19.30/95b8242ad729449fdb86591b9643c11463b61fa2/${sparc_tpl_mpi_ext}
export SGM_ROOT=${sparc_tpl_base}/sgm-19.30/00000000/${sparc_tpl_mpi_ext}
export SUPERLUDIST_ROOT=${sparc_tpl_base}/superlu_dist-5.4.0/a3121eaff44f7bf7d44e625c3b3d2a9911e58876/${sparc_tpl_mpi_ext}

export ATDM_CONFIG_USE_HWLOC=OFF
export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"
export ATDM_CONFIG_NETCDF_LIBS="-L${BOOST_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${PNETCDF_ROOT}/lib;-L${HDF5_ROOT}/lib;${BOOST_ROOT}/lib/libboost_program_options.a;${BOOST_ROOT}/lib/libboost_system.a;${NETCDF_ROOT}/lib/libnetcdf.a;${NETCDF_ROOT}/lib/libpnetcdf.a;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl;-lcurl"
export ATDM_CONFIG_SUPERLUDIST_INCLUDE_DIRS=${SUPERLUDIST_ROOT}/include
export ATDM_CONFIG_SUPERLUDIST_LIBS=${SUPERLUDIST_ROOT}/lib64/libsuperlu_dist.a
export ATDM_CONFIG_BINUTILS_LIBS="/usr/lib64/libbfd.so;/usr/lib64/libiberty.a"

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
export FC=mpif77
export F90=mpif90

export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE


#
# Run a script on the compute node send STDOUT and STDERR output to a file
# while also echo output to the console.  The primary purpose is to run the
# tests on the compute node.
#
# Usage:
#
#   atdm_run_script_on_comput_node <script_to_run> <output_file> \
#     [<timeout>] [<account>]
#
# If <timeout> and/or <account> are not given, then defaults are provided that
# work for the Jenkins driver process.
#
# In this case, sbatch is used to run the script but it also sends ouptut to
# STDOUT in real-time while it is running in addition to writing to the
# <outout_file>.  The job name for the sbatch script is taken from the env var
# 'ATDM_CONFIG_BUILD_NAME'.  This works for local builds since ATDM_CONFIG_BUILD_NAME.
#
# Note that you can pass in the script to run with arguments such as with
# "<some-script> <arg1> <arg2>" and it will work.  But note that this has to
# be bash script that 'sbatch' can copy and run form a temp location and it
# still has to work.  So the script has to use absolute directory paths, not
# relative paths or asume sym links, etc.
#
function atdm_run_script_on_compute_node {

  set +x

  script_to_run=$1
  output_file=$2
  timeout_input=$3
  account_input=$4

  echo
  echo "***"
  echo "*** atdm_run_script_on_compute_node '${script_to_run}' '${output_file}' '${timeout_input}' '${account_input}'"
  echo "***"
  echo

  if [ "${timeout_input}" == "" ] ; then
    timeout=1:30:00
  else
    timeout=${timeout_input}
  fi

  if [ "${account_input}" == "" ] ; then
    account=fy150090
  else
    account=${account_input}
  fi

  if [ -e $output_file ] ; then
    echo "Remove existing file $output_file"
    rm $output_file
  fi
  echo "Create empty file $output_file"
  touch $output_file

  echo
  echo "Running '$script_to_run' using sbatch in the background ..."
  set -x
  sbatch --output=$output_file --wait -N1 --time=${timeout} \
    -J $ATDM_CONFIG_BUILD_NAME --account=${account} ${script_to_run} &
  SBATCH_PID=$!
  set +x

  echo
  echo "Tailing output file $output_file in the background ..."
  set -x
  tail -f $output_file &
  TAIL_BID=$!
  set +x

  echo
  echo "Waiting for SBATCH_PID=$SBATCH_PID ..."
  wait $SBATCH_PID

  echo
  echo "Kill TAIL_BID=$TAIL_BID"
  kill -s 9 $TAIL_BID

  echo
  echo "Finished running ${script_to_run}!"
  echo

}

export -f atdm_run_script_on_compute_node

# NOTE: The above function is implemented in this way using 'sbatch' so that
# we can avoid using 'salloc' which is belived to cause ORTE errors.  But we
# still want to see live ouput from the script so that we can report it on
# Jenkins.  Therefore, the above approach is to use 'sbatch' and write its
# output to a known file-name.  Then, we use `tail -f` to print that file as
# it gets filled in from the 'sbatch' command.  The 'sbatch' command is run
# with --wait but is backgrouned to allow this to happen.  Then we wait for
# the 'sbatch' command to complete and then we kill the 'tail -f' command.
# That might seem overly complex but that gets the job done.
