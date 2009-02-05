#!/usr/bin/env sh
#

# This script changes file names for refactoring for explicit instantiation.
#
# You can run this multiple times on any base directory except the base
# directories for thyra mentioned below.
#
# You don't need to do any setup in order to run this script.
#
# WARNING!  Do not run this again from within the base thyra directories since
# it will mess things up.

# Get the directory for this scirpt which will give us the Trilinos base
# directory
#echo "0 = $0"
_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*\.sh/\1/g"`
#echo "_SCRIPT_DIR = $_SCRIPT_DIR"
_TRILINOS_HOME=$_SCRIPT_DIR/../../..
#echo "_TRILINOS_HOME = $_TRILINOS_HOME"

$_TRILINOS_HOME/commonTools/refactoring/token-replace-list-r \
  $_TRILINOS_HOME/packages/thyra/refactoring/explicit-instant.200902.token-list
