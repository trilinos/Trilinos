#!/usr/bin/env sh
#
# This script changes the names of the sub-vector and sub-multi-vector
# view classes.
#
# To run this script, you must have the variable $TRILINOS_HOME set
# to the base Trilinos directory so that
#
#    $TRILINOS_HOME/packages/rtop/refactoring
#
# is where this file is.  Trilinos/commonTools/refactoring must also
# be added to your path.
#

token-replace-list-r $TRILINOS_HOME/packages/rtop/refactoring/change-sub-multi-vector-view-names.20060311.token.list
