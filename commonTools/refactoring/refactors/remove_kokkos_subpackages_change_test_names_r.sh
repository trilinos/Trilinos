#!/bin/sh

# Script to update the names of Kokkos tests after the refactoring to remove
# subpackages.
#
# Run this script from the base directory and it will update all files.
#
# Warning!  Do not run this in a directory where the file
# remove_kokkos_subpackages_change_test_names.token-list will be found or this
# will change the names there too and destroy this script!
#
# Also, do not create a symbolic link to this script, just use a
# absolute or relative path from where it is defined in the
# Trilinos source tree.

# Get the directory for this scirpt which will give us the Trilinos base
# directory
_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*\.sh/\1/g"`
_TRILINOS_HOME=$_SCRIPT_DIR/../../..

# Run the replacements on all of the files found in subdirectories
find . -type f \
  -exec $_TRILINOS_HOME/commonTools/refactoring/string-replace-list.pl \
  $_SCRIPT_DIR/remove_kokkos_subpackages_change_test_names_r.token-list \
  '{}' '{}' ';'
