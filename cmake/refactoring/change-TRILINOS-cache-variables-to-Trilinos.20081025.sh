#!/usr/bin/env sh
#
# This script changes the names of some user cache varaibles from 'TRILINOS_' to
# 'Trilinos_' to be more consistent with other variable names.
#
# This script should be safe to run on code multiple times with no side effects

# Get the directory for this scirpt which will give us the Trilinos base
# directory
_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*\.sh/\1/g"`
_TRILINOS_HOME=$_SCRIPT_DIR/../..

find . -name CMakeLists.txt -exec $_TRILINOS_HOME/commonTools/refactoring/token-replace-list.pl \
  $_SCRIPT_DIR/change-TRILINOS-cache-variables-to-Trilinos.20081025.token.list {} {} \;

find . -name "*cmake" -exec $_TRILINOS_HOME/commonTools/refactoring/token-replace-list.pl \
  $_SCRIPT_DIR/change-TRILINOS-cache-variables-to-Trilinos.20081025.token.list {} {} \;
