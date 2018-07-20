################################################################################
#
# Set up env on mutrino for ATMD builds of Trilinos
#
# This source script gets the settings from the JOB_NAME var.
#
################################################################################

if [ "$ATDM_CONFIG_COMPILER" == "DEFAULT" ] ; then
  export ATDM_CONFIG_COMPILER=INTEL
fi

echo "Using mutrino compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

export ATDM_CONFIG_USE_NINJA=OFF
export ATDM_CONFIG_BUILD_COUNT=8
export OMP_NUM_THREADS=2

# Don't purge as this removes the Cray default environment
#module purge
module use /projects/EMPIRE/mutrino/tpls/hsw/modulefiles

# Use srun to run mpi jobs
export ATDM_CONFIG_MPI_EXEC="/opt/slurm/bin/srun"

# srun does not accept "-np" for # of processors
export ATDM_CONFIG_MPI_EXEC_NUMPROCS_FLAG="--ntasks"
export ATDM_CONFIG_MPI_PRE_FLAGS="--mpi=pmi2;--ntasks-per-node;36"
export ATDM_CONFIG_KOKKOS_ARCH=HSW
export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=16

if [ "$ATDM_CONFIG_COMPILER" == "INTEL" ]; then
    module load devpack/20180124/cray/7.6.2/intel/17.0.4
    module load gcc/4.9.3
    module load cmake/3.9.0
    export MPICXX=`which CC`
    export MPICC=`which cc`
    export MPIF90=`which ftn`

#    # Cray provides differently named wrappers
#    export CXX=`which CC`
#    export CC=`which cc`
    export ATDM_CONFIG_LAPACK_LIB="-mkl"
    export ATDM_CONFIG_BLAS_LIB="-mkl"
else
    echo "***"
    echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not supported on this system!"
    echo "***"
    return
fi

export ATDM_CONFIG_USE_HWLOC=OFF

export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"
export ATDM_CONFIG_NETCDF_LIBS="-L${BOOST_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${PNETCDF_ROOT}/lib;-L${HDF5_ROOT}/lib;${BOOST_ROOT}/lib/libboost_program_options.a;${BOOST_ROOT}/lib/libboost_system.a;${NETCDF_ROOT}/lib/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl;-lm"

export ATDM_CONFIG_MPI_POST_FLAG="-c 4"
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
# 'JOB_NAME'.  This works for local builds since JOB_NAME.
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
#  sbatch --output=$output_file --wait -N1 --time=${timeout} \
#    -J $JOB_NAME --account=${account} ${script_to_run} &
#  SBATCH_PID=$!
  sbatch --output=$output_file --wait -N1 --time=${timeout} \
    -J $JOB_NAME ${script_to_run} &
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

# NOTE: The above function is implemented in this way using 'sbatch' so that
# we can avoid using 'salloc' which is belived to cause ORTE errors.  But we
# still want to see live ouput from the script so that we can report it on
# Jenkins.  Therefore, the above approach is to use 'sbatch' and write its
# output to a known file-name.  Then, we use `tail -f` to print that file as
# it gets filled in from the 'sbatch' command.  The 'sbatch' command is run
# with --wait but is backgrouned to allow this to happen.  Then we wait for
# the 'sbatch' command to complete and then we kill the 'tail -f' command.
# That might seem overly complex but that gets the job done.
