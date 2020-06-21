function atdm_run_script_on_compute_node {
  script_to_run=$1
  output_file=$2
  ${script_to_run} 2>&1 | tee ${output_file}
}

export -f atdm_run_script_on_compute_node
