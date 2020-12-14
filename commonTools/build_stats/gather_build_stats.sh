#!/bin/bash
output_file=$1 ; shift
if [[ "${output_file}" == "" ]] ; then
  output_file=build_stats.csv
fi
find -name "*timing" -print -quit | xargs head -n1 > "${output_file}"
find -name "*timing" -print0 | xargs -0 tail -q -n1 >> "${output_file}"
