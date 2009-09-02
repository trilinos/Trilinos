#!/usr/bin/env sh
#
# This script changes a bunch of functions that started with a upper-case
# letter to start with a lower case letter.  This is to make it consistent
# with the Thyra naming convention in the "Thyra Coding and Doxygen
# Documentation Guildlines" document.
#
# Warning!  Do not run this in a directory where the
# lowercase-function-names.2007050515.token.list will be found or this will
# change the names there too!
#
# Warning!  This script changes the names of short functions so don't run this
# script on a large body of other code as it might result in changes that you
# don't want that are not related to this refactoring of rythmos.

# Get the directory for this scirpt which will give us the Trilinos base
# directory
_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*\.sh/\1/g"`
_TRILINOS_HOME=$_SCRIPT_DIR/../../..

# Run the replacements on all of the files found in subdirectories
$_TRILINOS_HOME/commonTools/refactoring/token-replace-list-r \
  $_SCRIPT_DIR/lowercase-function-names.2007050515.token.list
