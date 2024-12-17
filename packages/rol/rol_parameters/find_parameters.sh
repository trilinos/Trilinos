#!/bin/bash
# File: find_parameters.sh 

# Strip C/C++ style comments from source code
remove_comments() {
  perl -0777 -pe 's{//.*$}{}gm; s{/\*.*?\*/}{}gs'
}

# Check if a directory is provided as an argument, otherwise use the current directory
search_dir="${1:-.}"

# Verify that the directory exists
if [ ! -d "${search_dir}" ]; then
  echo "Error: Directory '${search_dir}' does not exist."
  echo "Usage: ${0} [directory_to_search]"
  exit 1
fi

# Updated pattern definitions
getset_pattern="(\.|\->)\s*((s|g)et|sublist)\s*\(\s*\""
sublist_pattern="(\.|\->)\s*sublist\s*\(\s*\""
pattern="${getset_pattern}|${sublist_pattern}"

# Function to generate find's -not -path arguments for excluded directories
generate_exclude_args() {
  local IFS=','
  local exclude_args=""
  for dir in $1; do
    exclude_args="${exclude_args} -not -path */${dir}/*"
  done
  echo ${exclude_args}
}

excluded_dirs="compatability,dynamic,interiorpoint,oed,sol,zoo"
exclude_args=$(generate_exclude_args "${excluded_dirs}")

find "${search_dir}" -type f -name '*.hpp' ${exclude_args} | while read -r file; do
  result=$(remove_comments < "${file}" | grep -E "${pattern}" || true)
  if [ -n "${result}" ]; then
    echo "$result" | while IFS= read -r line; do
      echo "${file}:${line}"
    done
  fi
done
