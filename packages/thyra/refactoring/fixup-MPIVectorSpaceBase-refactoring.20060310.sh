#!/usr/bin/env sh
#
# This script is just used to fixup a manual refactoring that I did
# by hand in breaking out the interface MPIVectorSpaceBaseBase.
#
# To run this script, you must have the variable $TRILINOS_HOME set
# to the base Trilinos directory so that
#
#    $TRILINOS_HOME/packages/thyra/refactoring
#
# is where this file is.  Trilinos/commonTools/refactoring must also
# be added to your path.
#

#token-replace-list-r $TRILINOS_HOME/packages/thyra/refactoring/fixup-MPIVectorSpaceBase-refactoring.20060310.token.list

# This is disabled since it really should not be used any more in any case.
