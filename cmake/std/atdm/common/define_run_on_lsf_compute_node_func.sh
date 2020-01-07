#TODO: s/SBATCH/LSF
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

#  if [[ "${timeout_input}" != "" ]] ; then
#    timeout=${timeout_input}
#  elif [[ "${ATDM_CONFIG_SBATCH_DEFAULT_TIMEOUT}" != "" ]] ; then
#    timeout=${ATDM_CONFIG_SBATCH_DEFAULT_TIMEOUT}
#  else
    timeout=4:OO
#  fi

  if [ "${account_input}" != "" ] ; then
    account=${account_input}
  elif [[ "${ATDM_CONFIG_SBATCH_DEFAULT_ACCOUNT}" != "" ]] ; then
    timeout=${ATDM_CONFIG_SBATCH_DEFAULT_ACCOUNT}
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
  #echo "bsub -J ${JOBNAME}-Test1 -W 06:00 -Is -n 16 -q ${QUEUE} ./test_submitted_command &> test.output" &> test_command
  bsub -J $ATDM_CONFIG_BUILD_NAME -W ${timeout} -Is ./${script_to_run} &> $output_file &
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
