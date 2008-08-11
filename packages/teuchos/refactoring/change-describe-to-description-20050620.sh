#!/bin/sh

# This script can be used to automatically change the name of the function
# describe() to description() (in Teuchos::Describable and its subclasses).
# To use this script you must set the environment
# variable TRILINOS_HOME to the path of your Trilinos base directory and you
# must add $TRILINOS_HOME/commonTools/refactoring to your path.

# Run this script from the base directory of any code that you would like to
# upgrade.  This scrpit can be run mutiple times with no side effects.

string-replace-list-r $TRILINOS_HOME/packages/teuchos/refactoring/change-describe-to-description-20050620.string-list

