#!/bin/sh

# This script can be used to automatially update code to change
# from Teuchos::RefCountPtr to Teuchos::RCP.
#
# Run this script from any base directory where you want to update
# the files.  Note, this will update all files, except those
# in a directory CVS.
#
# Warning!  Do not run this in a directory where the
# file change-RefCountPtr-to-RCP-20070619.token-list will be found or
# this will change the names there too and destroy this script!
#
# Also, do not create a symbolic link to this script, just use a
# absolute or relative path from where it is defined in the
# Trilinos source tree.

# Get the directory for this scirpt which will give us the Trilinos base
# directory
_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*\.sh/\1/g"`
_TRILINOS_HOME=$_SCRIPT_DIR/../../..

# Run the replacements on all of the files found in subdirectories
$_TRILINOS_HOME/commonTools/refactoring/token-replace-list-r \
  $_SCRIPT_DIR/change-RefCountPtr-to-RCP-20070619.token-list
