#
# Run a script on the compute node send STDOUT and STDERR output to a file
# while also echo output to the console.  The primary purpose is to run a
# arbitrary script on a single compute node.
#
# Usage:
#
#   atdm_run_script_on_compute_node <script_to_run> <output_file> [<timeout>]
#
# If <timeout> is not given, then it is set to the env var
# ATDM_CONFIG_BSUB_DEFAULT_TIMEOUT if != "".  Otherwise, a default is set
# internally.
#
# In this case, "lalloc 1 <script_tot_run>" is used to run the script but it
# also sends ouptut to STDOUT in real-time while it is running in addition to
# writing to the <outout_file>.  The SLURM job name for the sbatch script is
# taken from the env var 'ATDM_CONFIG_BUILD_NAME'.  This works for local
# builds as well since ATDM_CONFIG_BUILD_NAME is always set by the
# atdm/load-env.sh script.
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

  echo
  echo "***"
  echo "*** atdm_run_script_on_compute_node '${script_to_run}' '${output_file}' '${timeout_input}''"
  echo "***"
  echo

  if [[ "${timeout_input}" != "" ]] ; then
    timeout=${timeout_input}
  elif [[ "${ATDM_CONFIG_BSUB_DEFAULT_TIMEOUT}" != "" ]] ; then
    timeout=${ATDM_CONFIG_BSUB_DEFAULT_TIMEOUT}
  else
    timeout=04:00
  fi

  if [ -e $output_file ] ; then
    echo "Remove existing file $output_file"
    rm $output_file
  fi
  echo "Create empty file $output_file"
  touch $output_file

  echo
  echo "Running '$script_to_run' on compute node in the background ..."
  set -x
  lalloc 1 -W ${timeout} ${script_to_run} &> $output_file &
  BSUB_PID=$!
  set +x

  echo
  echo "Tailing output file $output_file in the background ..."
  set -x
  tail -f $output_file &
  TAIL_PID=$!
  set +x

  echo
  echo "Waiting for BSUB_PID=$BSUB_PID ..."
  wait $BSUB_PID

  echo
  echo "Kill TAIL_BID=$TAIL_BID"
  kill -s 9 $TAIL_PID

  echo
  echo "Finished running ${script_to_run}!"
  echo

}

export -f atdm_run_script_on_compute_node

# NOTE: The above function is implemented in this way using 'BSUB' to improve
# test run times (submitting a job for each test is slow). We  still want to
# see live ouput from the script so that we can report it on
# Jenkins and see it running in real-time on Jenins.  Therefore, the above
# approach is to use 'BSUB' and write its output to a known file-name.
# Then, we use `tail -f` to print that file as it gets filled in from the
# 'BSUB' command.  We wait for the 'BSUB' command  to complete and then we kill
# the 'tail -f' command.  That might seem overly  complex but that gets the job
# done.
