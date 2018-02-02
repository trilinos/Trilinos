#!/bin/bash -e

#
# This script can be run from any directory as:
#
#   <some-base-dir>/Trilinos/cmake/update-tribits-snapshot.sh [other options]
#
# to update the snapshot of Trilinos/TriBITS/tribits/ into
# Trilios/cmake/tribits/.  This snapshot script should only be run whil on the
# 'tribits_github_snapshot' branch of Trilinos.  And then that branch is
# merged back to the main Trilinos 'develop' branch.  Details on how to set up
# for and execute snapshotting TriBITS are given in:
#
#  http://trac.trilinos.org/wiki/TriBITSTrilinosDev
#
# NOTE: Do not attempt to update the snapshot of TriBITS without consulting
# that documentation!
#

_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*\.sh/\1/g"`
echo "_SCRIPT_DIR = '${_SCRIPT_DIR}'"
_TRILINOS_HOME=$_SCRIPT_DIR/..

# Update the snapshot
cd $_TRILINOS_HOME/cmake/tribits/
../../TriBITS/tribits/snapshot_tribits.py --clean-ignored-files-orig-dir "$@"
