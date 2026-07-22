#!/bin/sh
#
# Refactor a TriBITS project's CMake files for the change from
# find_package(PythonInterp) to find_package(Python3) (see CHANGELOG.md entry
# for 2024-10-08).
#
# Usage:
#
#   to-python3.sh <base-dir>
#
# NOTES:
# * This is safe to run multiple times on the same set of files.
# * This is safe to run from the base project source tree.
# * This will exclude any files under '.git/'
# * This will exclude any files under 'tribits/' directories (so it will not mess up TriBITS)
#

_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/to-python3[.]sh/\1/g"`
baseDir=$1
find ${baseDir} -type f \
  \( -name CMakeLists.txt -or -name "*.cmake" -or -name "*.cmake.in" -or -name "*.rst" \) \
  ! -wholename "*/.git/*" ! -wholename  "*/tribits/*" \
  -exec $_SCRIPT_DIR/token-replace.py -t PYTHON_EXECUTABLE -r Python3_EXECUTABLE -f {} ';'
