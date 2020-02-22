if [[ `type -t atdm_run_script_on_compute_node`"" != 'function' ]]; then

  # Create a dummy function that just runs the script, writes STDOUT to a file
  # and sends STDOUT to STDOUT.

  function atdm_run_script_on_compute_node {
    script_to_run=$1
    output_file=$2
    ${script_to_run} 2>&1 | tee ${output_file}
  }

  export -f atdm_run_script_on_compute_node

else
  echo "Warning: atdm_run_script_on_compute_node is already set!"
fi
