################################################################################
#
# Set up env on mutrino for ATMD builds of Trilinos
#
# This source script gets the settings from the ATDM_CONFIG_BUILD_NAME var.
#
################################################################################

if [ "$ATDM_CONFIG_COMPILER" == "DEFAULT" ] ; then
  export ATDM_CONFIG_COMPILER=INTEL
fi

if [ "$ATDM_CONFIG_KOKKOS_ARCH" == "DEFAULT" ] ; then
  export ATDM_CONFIG_KOKKOS_ARCH=HSW
fi

if [ "$ATDM_CONFIG_KOKKOS_ARCH" != "HSW" ] && [ "$ATDM_CONFIG_KOKKOS_ARCH" != "KNL" ]; then
  echo
  echo "***"
  echo "*** ERROR: KOKKOS_ARCH=$ATDM_CONFIG_KOKKOS_ARCH is not a valid option on this system."
  echo "*** '$ATDM_CONFIG_KOKKOS_ARCH' appears in $ATDM_CONFIG_BUILD_NAME which then sets the KOKKOS_ARCH."
  echo "*** On Mutrino 'HSW' and 'KNL' are the only valid KOKKOS_ARCH settings."
  echo "*** If no KOKKOS_ARCH is specified then 'HSW' will be used by default."
  echo "***"
  return
fi

echo "Using mutrino compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE and KOKKOS_ARCH=$ATDM_CONFIG_KOKKOS_ARCH"

export ATDM_CONFIG_USE_NINJA=ON
export ATDM_CONFIG_BUILD_COUNT=8
# NOTE: We are building on the one 'mutrino' login node so we can't use a lot
# of cores to build.

# Use srun to run mpi jobs
export ATDM_CONFIG_MPI_EXEC="/opt/slurm/bin/srun"

# srun does not accept "-np" for # of processors
export ATDM_CONFIG_MPI_EXEC_NUMPROCS_FLAG="--ntasks"
export ATDM_CONFIG_MPI_PRE_FLAGS="--mpi=pmi2;--ntasks-per-node;36"

export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=8
# NOTE: Using -j8 instead of -j16 for ctest is to try to avoid 'srun' "Job
# <jobid> step creation temporarily disabled" failures on 'mutrino' (see
# TRIL-214).

if [ "$ATDM_CONFIG_COMPILER" == "INTEL" ] && [ "$ATDM_CONFIG_KOKKOS_ARCH" == "HSW"  ]; then
    module use /projects/EMPIRE/mutrino/tpls/hsw/modulefiles
    export OMP_NUM_THREADS=2
    export ATDM_CONFIG_MPI_POST_FLAGS="-c 4"
elif [ "$ATDM_CONFIG_COMPILER" == "INTEL" ] && [ "$ATDM_CONFIG_KOKKOS_ARCH" == "KNL"  ]; then
    module use /projects/EMPIRE/mutrino/tpls/knl/modulefiles
    export SLURM_TASKS_PER_NODE=16
    export OMP_NUM_THREADS=2
    export OMP_PLACES=threads
    export OMP_PROC_BIND=spread
    export ATDM_CONFIG_MPI_POST_FLAGS="--hint=nomultithread;-c 4"
    export ATDM_CONFIG_SBATCH_OPTIONS="-p knl -C cache --hint=multithread"
else
    echo
    echo "***"
    echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER and KOKKOS_ARCH=$ATDM_CONFIG_KOKKOS_ARCH is not"
    echo "*** a supported combination on this system!"
    echo "*** Combinations that are supported: "
    echo "*** > Intel compiler with KOKKOS_ARCH=HSW"
    echo "*** > Intel compiler with KOKKOS_ARCH=KNL"
   echo "***"
    return
fi

# Load the modules (can't purge)
module load devpack/20180124/cray/7.6.2/intel/17.0.4
module load gcc/4.9.3
module load cmake/3.9.0

# No RPATH for static builds
export ATDM_CONFIG_CMAKE_SKIP_INSTALL_RPATH=ON
# ToDo: Make above contingent on 'static' or 'shared' 

# Use manually installed cmake and ninja to allow usage of ninja and
# all-at-once mode
export PATH=/projects/netpub/atdm/cmake-3.11.4/bin:/projects/netpub/atdm/ninja-1.8.2/bin:$PATH

# Set MPI wrappers
export MPICXX=`which CC`
export MPICC=`which cc`
export MPIF90=`which ftn`

# Cray provides differently named wrappers
export ATDM_CONFIG_LAPACK_LIBS="-mkl"
export ATDM_CONFIG_BLAS_LIBS="-mkl"

export ATDM_CONFIG_USE_HWLOC=OFF

export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"
export ATDM_CONFIG_NETCDF_LIBS="-L${BOOST_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${PNETCDF_ROOT}/lib;-L${HDF5_ROOT}/lib;${BOOST_ROOT}/lib/libboost_program_options.a;${BOOST_ROOT}/lib/libboost_system.a;${NETCDF_ROOT}/lib/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl;-lm"

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

  if [ "${output_file}" == "" ] ; then
    echo "ERROR: must specify output file as second argument!"
    return
  fi

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

  if [ "$ATDM_CONFIG_KOKKOS_ARCH" == "KNL" ] ; then
    ATDM_CONFIG_SBATCH_EXTRA_ARGS="-p knl -C cache --hint=multithread"
  elif [ "$ATDM_CONFIG_KOKKOS_ARCH" == "HSW" ] ; then
    ATDM_CONFIG_SBATCH_EXTRA_ARGS=
  else
   echo
   echo "***"
   echo "*** ERROR: Invalid value ATDM_CONFIG_KOKKOS_ARCH=${ATDM_CONFIG_KOKKOS_ARCH} specified!" 
   echo "***"
   return
  fi
 
  echo
  echo "Running '$script_to_run' using sbatch in the background ..."
  set -x
  sbatch --output=$output_file --wait -N1 ${ATDM_CONFIG_SBATCH_EXTRA_ARGS} \
    --time=${timeout} -J $ATDM_CONFIG_BUILD_NAME ${script_to_run} &
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
