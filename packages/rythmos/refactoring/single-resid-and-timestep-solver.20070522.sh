#!/usr/bin/env sh
#
# This script updates code for moving Thyra::SingleResidSSDAEModelEvaluator
# to Rytymos::SingleResidualModelEvaluator and for moving
# Thyra::TimeStepNewtonNonlinearSolver to Rythmos::TimeStepNonlinearSolver.
#
# Warning!  Do not run this in a directory where the
# single-resid-and-timestep-solver.20070522.token.list will be found or this will
# change the names there too!
#
# This script should be safe to run on code multiple times with no side effects

# Get the directory for this scirpt which will give us the Trilinos base
# directory
_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*\.sh/\1/g"`
_TRILINOS_HOME=$_SCRIPT_DIR/../../..

# Run the replacements on all of the files found in subdirectories
$_TRILINOS_HOME/commonTools/refactoring/token-replace-list-r \
  $_SCRIPT_DIR/single-resid-and-timestep-solver.20070522.token.list
