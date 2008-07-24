#!/usr/bin/env sh
#

# This script changes the name of the class ReductTargetScalar to
# DefaultReductTarget.  It does not obey namespace names so be a little
# careful.

# You can run this script on the same source code as much as you want but just
# be careful not to run it in this directory or it will change the names of
# the token replace list.

# Get the directory for this scirpt which will give us the Trilinos base
# directory
_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*\.sh/\1/g"`
_TRILINOS_HOME=$_SCRIPT_DIR/../../..

$_TRILINOS_HOME/commonTools/refactoring/token-replace-list-r \
  $_TRILINOS_HOME/packages/rtop/refactoring/change-ReductTargetScalar-to-DefaultReductTarget.20080721.token.list
