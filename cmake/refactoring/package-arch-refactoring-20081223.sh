#!/usr/bin/env sh
#

# This script changes the names of include CMake files, macro and
# function names, etc.

#
# This script should be safe to run on code multiple times with no side effects

# Get the directory for this scirpt which will give us the Trilinos base
# directory
_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*\.sh/\1/g"`
_TRILINOS_HOME=$_SCRIPT_DIR/../..

find . -name CMakeLists.txt -exec $_TRILINOS_HOME/commonTools/refactoring/token-replace-list.pl \
  $_SCRIPT_DIR/package-arch-refactoring-20081223.token.list {} {} \;

find . -name "*cmake" -exec $_TRILINOS_HOME/commonTools/refactoring/token-replace-list.pl \
  $_SCRIPT_DIR/package-arch-refactoring-20081223.token.list {} {} \;
