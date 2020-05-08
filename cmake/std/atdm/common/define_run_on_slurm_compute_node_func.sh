#
# Run a script on the compute node send STDOUT and STDERR output to a file
# while also echo output to the console.  The primary purpose is to run the
# tests on the compute node.
#
# Usage:
#
#   atdm_run_script_on_compute_node <script_to_run> <output_file> \
#     [<timeout>] [<account>]
#
# If <timeout> is not given, then it is set to the env var
# ATDM_CONFIG_SBATCH_DEFAULT_TIMEOUT if != "".  Otherwise, a default is set
# internally.
#
# If <account> is not given, then it is set to the env var
# ATDM_CONFIG_SBATCH_DEFAULT_ACCOUNT if != "".  Otherwise, a default is set
# internally.
#
# The value of ATDM_CONFIG_SBATCH_EXTRA_ARGS is added to sbatch as follows:
#   sbatch --output=$output_file --wait -N1 \
#   ${ATDM_CONFIG_SBATCH_EXTRA_ARGS} --time=${timeout} \
#   -J $ATDM_CONFIG_BUILD_NAME --account=${account} ${script_to_run} &
#
# In this case, sbatch is used to run the script but it also sends ouptut to
# STDOUT in real-time while it is running in addition to writing to the
# <outout_file>.  The SLURM job name for the sbatch script is taken from the
# env var 'ATDM_CONFIG_BUILD_NAME'.  This works for local builds as well since
# ATDM_CONFIG_BUILD_NAME is always set by the atdm/load-env.sh script.
#
# Note that you can pass in the script to run with arguments such as with
# "<some-script> <arg1> <arg2>" and it will work.  But note that this script
# has to be bash script that 'sbatch' can copy and run from a temp location
# and it still has to work.  Therefore, this script has to use absolute
# directory paths, not relative paths or asume sym links, etc (using env vars
# to know what to run and where).
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

  if [[ "${timeout_input}" != "" ]] ; then
    timeout=${timeout_input}
  elif [[ "${ATDM_CONFIG_SBATCH_DEFAULT_TIMEOUT}" != "" ]] ; then
    timeout=${ATDM_CONFIG_SBATCH_DEFAULT_TIMEOUT}
  else
    timeout=1:30:00
  fi

  if [ "${account_input}" != "" ] ; then
    account=${account_input}
  elif [[ "${ATDM_CONFIG_SBATCH_DEFAULT_ACCOUNT}" != "" ]] ; then
    account=${ATDM_CONFIG_SBATCH_DEFAULT_ACCOUNT}
  else
    account=fy150090
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
  sbatch --output=$output_file --wait -N1 ${ATDM_CONFIG_SBATCH_EXTRA_ARGS} \
    --time=${timeout} -J $ATDM_CONFIG_BUILD_NAME --account=${account} ${script_to_run} &
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
# Jenkins and see it running in real-time on Jenins.  Therefore, the above
# approach is to use 'sbatch' and write its output to a known file-name.
# Then, we use `tail -f` to print that file as it gets filled in from the
# 'sbatch' command.  The 'sbatch' command is run with --wait but is
# backgrouned to allow this to happen.  Then we wait for the 'sbatch' command
# to complete and then we kill the 'tail -f' command.  That might seem overly
# complex but that gets the job done.
