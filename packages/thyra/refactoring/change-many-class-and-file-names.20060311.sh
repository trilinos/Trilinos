#!/usr/bin/env sh
#
# This script changes the names of many different classes and therefore
# also file names.
#
# To run this script, you must have the variable $TRILINOS_HOME set
# to the base Trilinos directory so that
#
#    $TRILINOS_HOME/packages/thyra/refactoring
#
# is where this file is.  Trilinos/commonTools/refactoring must also
# be added to your path.
#

token-replace-list-r $TRILINOS_HOME/packages/thyra/refactoring/change-many-class-and-file-names.20060311.token.list

