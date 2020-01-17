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
  echo "Running '$script_to_run' using BSUB in the background ..."
  set -x
  bsub -J $ATDM_CONFIG_BUILD_NAME -W ${timeout} -Is ${script_to_run} &> $output_file &
  BSUB_PID=$!
  set +x

  echo
  echo "Tailing output file $output_file in the background ..."
  set -x
  tail -f $output_file &
  TAIL_BID=$!
  set +x

  echo
  echo "Waiting for BSUB_PID=$BSUB_PID ..."
  wait $BSUB_PID

  echo
  echo "Kill TAIL_BID=$TAIL_BID"
  kill -s 9 $TAIL_BID

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
