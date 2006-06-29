#!/usr/bin/env sh
#
# This script changes the names of classes and file names from
# for the replacement of the separate Serial and MPI support software
# to a single set of Spmd classes based on the new Tuechos::Comm interface.
#
# To run this script, you must have the variable $TRILINOS_HOME set
# to the base Trilinos directory so that
#
#    $TRILINOS_HOME/packages/thyra/refactoring
#
# is where this file is.  Trilinos/commonTools/refactoring must also
# be added to your path.
#

token-replace-list-r $TRILINOS_HOME/packages/thyra/refactoring/spmd-support-refactoring.20060627.token.list
