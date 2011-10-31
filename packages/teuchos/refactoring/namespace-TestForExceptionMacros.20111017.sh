#!/bin/sh

# This script is run to update all client code to replace the non-namespaced
# macros in Teuchos_TestForException.hpp to the namespaced macros in
# Teuchos_Assert.hpp.
#
# NOTE: This script is safe to run multiple times on the same directory and is
# even safe to run on the teuchos directory that contain these files (because
# of the file ignore list).

# Get the directory for this script which will give us the Trilinos base
# directory
_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*\.sh/\1/g"`
_TRILINOS_HOME=$_SCRIPT_DIR/../../..

# Run the replacements on all of the files found in subdirectories
$_TRILINOS_HOME/commonTools/refactoring/token-replace-list-r \
  $_SCRIPT_DIR/namespace-TestForExceptionMacros.20111017.token-list \
  0 \
  $_SCRIPT_DIR/namespace-TestForExceptionMacros.20111017.ignore-files-list
