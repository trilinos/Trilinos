#!/usr/bin/env sh
#
# This script changes from ETransp to EOpTransp to disambiguate with
# Teuchos::ETransp since it caused trouble with Doxygen (see bug 3902)
#
# You can run this multiple times on any base directory except the base
# directories for thyra mentioned below.
#
# You don't need to do any setup in order to run this script.
#
# WARNING!  Do not run this again from within the base thyra or thyra/src
# directories since it will mess things up.

# Get the directory for this scirpt which will give us the Trilinos base
# directory
_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*\.sh/\1/g"`
_TRILINOS_HOME=$_SCRIPT_DIR/../../..

$TRILINOS_HOME/commonTools/refactoring/token-replace-list-r \
  $TRILINOS_HOME/packages/thyra/refactoring/change-Thyra-ETransp-to-EOpTransp.20080520.token_list
