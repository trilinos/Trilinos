#!/bin/sh

# This script can be used to automatically change the name of the macro
# _DEBUG (which is used by other software with bad sideeffects) to
# TEUCHOS_DEBUG.   To use this script you must set the environment
# variable TRILINOS_HOME to the path of your Trilinos base directory and you
# must add $TRILINOS_HOME/commonTools/refactoring to your path.

# Run this script from the base directory of any code that you would like to
# upgrade.  This scrpit can be run mutiple times with no side effects.

token-replace-list-r $TRILINOS_HOME/packages/teuchos/refactoring/change-_DEBUG-to-TEUCHOS_DEBUG-20060620.token-list
