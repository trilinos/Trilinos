#!/bin/sh

# This script changes the include names of several files to standarize on the
# Tribits* prefix.
#
# NOTE: This script is safe to run multiple times on the same directory and is
# even safe to run on the teuchos directory that contain these files (because
# of the file ignore list).

# Get the directory for this script which will give us the Trilinos base
# directory
_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*\.sh/\1/g"`
_TRILINOS_HOME=$_SCRIPT_DIR/../..

find . -name "CMakeLists.txt" \
  -exec $_TRILINOS_HOME/commonTools/refactoring/token-replace-list.pl \
    $_SCRIPT_DIR/TrilinosCreateClientTemplateHeaders-refactor.20111202.token-list \
   {} {} 0 \
   $_SCRIPT_DIR/TrilinosCreateClientTemplateHeaders-refactor.20111202.ignore-files-list \
   \;

find . -name "*.cmake" \
  -exec $_TRILINOS_HOME/commonTools/refactoring/token-replace-list.pl \
    $_SCRIPT_DIR/TrilinosCreateClientTemplateHeaders-refactor.20111202.token-list \
   {} {} 0 \
   $_SCRIPT_DIR/TrilinosCreateClientTemplateHeaders-refactor.20111202.ignore-files-list \
   \;
