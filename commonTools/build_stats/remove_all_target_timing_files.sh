#!/bin/bash -e

base_dir=$1 ; shift
if [[ "${base_dir}" == "" ]] ; then
  "Error, must provide base dir as first argument!"
  exit 1
fi

find ${base_dir} -name "*.timing" -exec rm {} \;
#find ${base_dir} -name "*.timing" -exec ls -l {} \;
