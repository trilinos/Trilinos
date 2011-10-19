#!/bin/sh

# This script is run to update all client code to replace the non-namespaced
# macros in Teuchos_TestForException.hpp to the namespaced macros in
# Teuchos_Assert.hpp.
#
# WARNING: Do not run this script on the teuchos directory itself or it will
# mess up the core files Teuchos_Assert.*pp and Teuchos_TestForException.hpp.

# Get the directory for this scirpt which will give us the Trilinos base
# directory
_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*\.sh/\1/g"`
_TRILINOS_HOME=$_SCRIPT_DIR/../../..

# Run the replacements on all of the files found in subdirectories
$_TRILINOS_HOME/commonTools/refactoring/token-replace-list-r \
  $_SCRIPT_DIR/namespace-TestForExceptionMacros.20111017.token-list
