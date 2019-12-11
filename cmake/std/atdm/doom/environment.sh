################################################################################
#
# Set up env on the doom system for ATMD builds of Trilinos
#
# This source script gets the settings from the JOB_NAME var.
#
################################################################################

if   [[ "$ATDM_CONFIG_COMPILER" == "CUDA" ]] \
  || [[ "$ATDM_CONFIG_COMPILER" == "DEFAULT" ]] \
  ; then
  ATDM_CONFIG_COMPILER="CUDA-9.2_GNU-6.3.1_OPENMPI-2.1.1"
fi

#TODO: jfrye 
if [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "DEFAULT" ]] ; then
  ATDM_CONFIG_KOKKOS_ARCH="BDW,Pascal60"
else
  echo
  echo "***"
  echo "*** ERROR: Specifying KOKKOS_ARCH is not supported on doom builds"
  echo "*** remove '$ATDM_CONFIG_KOKKOS_ARCH' from JOB_NAME=$JOB_NAME"
  echo "***"
  return
fi

echo "Using doom compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

if [[ ! -d /tmp/${USER} ]]; then
  # Need to manually create this directory or CMake compiler check fails.
  # This must be a defect in nvcc.
  mkdir /tmp/${USER}
fi

export ATDM_CONFIG_ENABLE_SPARC_SETTINGS=ON
export ATDM_CONFIG_USE_NINJA=ON

# Get ATDM_CONFIG_NUM_CORES_ON_MACHINE for this machine
source $ATDM_SCRIPT_DIR/utils/get_num_cores_on_machine.sh

if [ "$ATDM_CONFIG_NUM_CORES_ON_MACHINE" -gt "16" ] ; then
  export ATDM_CONFIG_MAX_NUM_CORES_TO_USE=16
  # NOTE: We get links crashing if we try to use to many processes.  ToDo: We
  # should limit the number of processes that ninja uses to link instead of
  # reducing the overrall parallel build level like this.
else
  export ATDM_CONFIG_MAX_NUM_CORES_TO_USE=$ATDM_CONFIG_NUM_CORES_ON_MACHINE
fi

export ATDM_CONFIG_BUILD_COUNT=$ATDM_CONFIG_MAX_NUM_CORES_TO_USE
# NOTE: Use as many build processes and there are cores by default.

module purge

# For now, turn on warnings by default:
if [[ "${ATDM_CONFIG_ENABLE_STRONG_WARNINGS}" == "" ]] ; then
  export ATDM_CONFIG_ENABLE_STRONG_WARNINGS=1
fi

if [[ "$ATDM_CONFIG_COMPILER" == "CUDA-9.2_GNU-6.3.1_OPENMPI-2.1.1" ]] ; then
  module load sparc-dev/cuda-9.2.88_gcc-6.3.1_openmpi-2.1.1
  unset OMP_NUM_THREADS  # SPARC module sets these and we must unset!
  unset OMP_PROC_BIND
  unset OMP_PLACES
  export OMPI_CXX=`which g++`
  export OMPI_CC=`which gcc`
  export OMPI_FC=`which gfortran`
  export MPICC=`which mpicc`
  export MPICXX=`which mpicxx`
  export MPIF90=`which mpif90`
  export ATDM_CONFIG_MKL_ROOT=${CBLAS_ROOT}
  export ATDM_CONFIG_MPI_EXEC=mpirun

  export OMPI_CXX=${ATDM_CONFIG_NVCC_WRAPPER}
  if [ ! -x "$OMPI_CXX" ]; then
      echo "No nvcc_wrapper found"
      return
  fi
  # some Trilinos tests require this to run correctly
  export CUDA_LAUNCH_BLOCKING=1
  export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1
  export KOKKOS_NUM_DEVICES=2

else
  echo
  echo "***"
  echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not supported on this system!"
  echo "***"
fi

# ToDo: Update above to only load the compiler and MPI moudles and then
# directly set <TPL_NAME>_ROOT to point to the right TPLs.  This is needed to
# avoid having people depend on the SPARC modules.

# Set parallel test level based on OpenMP or not
if [[ "$ATDM_CONFIG_NODE_TYPE" == "OPENMP" ]] ; then
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=$(($ATDM_CONFIG_MAX_NUM_CORES_TO_USE/2))
  #export OMP_NUM_THREADS=2
else
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=$(($ATDM_CONFIG_MAX_NUM_CORES_TO_USE/2))
fi
# NOTE: Above, we use 1/2 as many parallel processes as cores on the machine
# to be safe.  Also, we need to set OMP_* env vars here because the SPARC
# modules change them!

export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=4

# Use updated Ninja (not provided by sparc-dev module above)
module load sems-env
module load atdm-env
module load atdm-ninja_fortran/1.7.2

export ATDM_CONFIG_USE_HWLOC=OFF

export ATDM_CONFIG_BINUTILS_LIBS="/usr/lib64/libbfd.so;/usr/lib64/libiberty.a"
# NOTE: Above, we have to explicitly set the libs to use libbdf.so instead of
# libbdf.a because the former works and the latter does not and TriBITS is set
# up to only find static libs by default!

# BLAS and LAPACK

atdm_config_add_libs_to_var ATDM_CONFIG_BLAS_LIBS ${ATDM_CONFIG_MKL_ROOT}/mkl/lib/intel64 .so \
  mkl_intel_lp64 mkl_intel_thread mkl_core

atdm_config_add_libs_to_var ATDM_CONFIG_BLAS_LIBS ${ATDM_CONFIG_MKL_ROOT}/../../lib/intel64 .so \
  iomp5

#atdm_config_add_libs_to_var ATDM_CONFIG_BLAS_LIBS ${ATDM_CONFIG_MKL_ROOT}/../compiler/lib/intel64_lin .so \
#  iomp5

export ATDM_CONFIG_LAPACK_LIBS=${ATDM_CONFIG_BLAS_LIBS}

# Boost

atdm_config_add_libs_to_var ATDM_CONFIG_BOOST_LIBS ${BOOST_ROOT}/lib .a \
  boost_program_options boost_system

# NOTE: Above, the SPARC-installed TPLs only have *.a files.  There are no
# *.so files.

# HDF5 and Netcdf

# NOTE: HDF5_ROOT and NETCDF_ROOT should already be set in env from above
# module loads!

# However, set the direct libs for HDF5 and NetCDF in case we use that option
# for building (see env var ATDM_CONFIG_USE_SPARC_TPL_FIND_SETTINGS).

export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"

export ATDM_CONFIG_NETCDF_LIBS="-L${BOOST_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${SEMS_PNETCDF_ROOT}/lib;-L${HDF5_ROOT}/lib;${BOOST_ROOT}/lib/libboost_program_options.a;${BOOST_ROOT}/lib/libboost_system.a;${NETCDF_ROOT}/lib/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl;-lcurl"

# SuperLUDist
if [[ "${ATDM_CONFIG_SUPERLUDIST_INCLUDE_DIRS}" == "" ]] ; then
  # Set the default which is correct for all of the new TPL builds
  export ATDM_CONFIG_SUPERLUDIST_INCLUDE_DIRS=${SUPERLUDIST_ROOT}/include
  export ATDM_CONFIG_SUPERLUDIST_LIBS=${SUPERLUDIST_ROOT}/lib64/libsuperlu_dist.a
fi


export ATDM_CONFIG_MPI_PRE_FLAGS="--bind-to;core:overload-allowed"

# Finished!
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
