#!/usr/bin/env sh
#
# This script performs any refactoring needed to update code for the reorganization
# of thyra into rtop and thyra.
#
# To run this script, you must have the variable $TRILINOS_HOME set
# to the base Trilinos directory so that
#
#    $TRILINOS_HOME/packages/thyra/refactoring
#
# is where this file is.  Trilinos/commonTools/refactoring must also
# be added to your path.
#

token-replace-list-r $TRILINOS_HOME/packages/thyra/refactoring/reorganize-thyra.20050201.token.list
