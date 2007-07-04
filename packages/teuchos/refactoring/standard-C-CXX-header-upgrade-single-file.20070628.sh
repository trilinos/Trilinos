#!/bin/sh

# This script can be used to automatically upgrade any C++ source file that
# uses C-style headers like <math.h> and expects standard library symbols like
# in the global namespace liek sqrt(...) and replace them with standard C++
# headers like <cmath> and puts in explicit namespace qualification
# std::sqrt(...).
#
# Run this script on any C++ header or source file that you need to upgrade as:
#
#  $TRILINOS_HOME/teuchos/refactoring/standard-C-CXX-header-upgrade-single-file.20070628.sh file
#
# Warning! Do not run this in on C files!
#
# Also, do not create a symbolic link to this script, just use a
# absolute or relative path from where it is defined in the
# Trilinos source tree.

# Get the directory for this scirpt which will give us the Trilinos base
# directory
_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*\.sh/\1/g"`
_TRILINOS_HOME=$_SCRIPT_DIR/../../..

$_TRILINOS_HOME/commonTools/refactoring/token-replace-list.pl \
  $_SCRIPT_DIR/standard-C-CXX-header-upgrade.20070628.token-list $1 $1 -v
$_TRILINOS_HOME/commonTools/refactoring/string-replace-list.pl \
  $_SCRIPT_DIR/standard-C-CXX-header-upgrade.20070628.string-list $1 $1
sed  -i 's/\([^a-zA-Z_]\)Teuchosstd::map\([^a-zA-Z0-9_]\)/\1Teuchos::map\2/g' $1

# Note: The bare token 'map' is not replaced with 'std::map' since the concept
# of a map is very common in Trilinos and this replacment causes more problems
# than it fixes!
