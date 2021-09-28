#!/bin/sh
#
# Usage:
#
#   lower-case-cmake-tree.sh <base-dir>
#
_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/lower-case-cmake-tree[.]sh/\1/g"`
baseDir=$1
find ${baseDir} -type f \
  \( -name CMakeLists.txt -or -name "*.cmake" -or -name "*.cmake.in" -or -name "*.rst" \) \
  -exec $_SCRIPT_DIR/../python_utils/lower_case_cmake.py {} ';'
