#!/usr/bin/env sh
#
# This script changes the names of some of the wrapper functions for
# the Epetra/Thyra adapters which followed the SPMD refactoring in the
# Thyra support software.
#
# To run this script, you must have the variable $TRILINOS_HOME set
# to the base Trilinos directory so that
#
#    $TRILINOS_HOME/packages/epetra/thyra/refactoring
#
# is where this file is.  Trilinos/commonTools/refactoring must also
# be added to your path.
#

token-replace-list-r $TRILINOS_HOME/packages/epetra/thyra/refactoring/spmd-support-refactoring.20060627.token.list
