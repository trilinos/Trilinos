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

# Get the merge commits being pulled in and amend the commit message
newTribitsMergeCommits=$(./cmake/get-tribits-mainline-merge-commits.sh)
#echo "newTribitsMergeCommits:
#
#${newTribitsMergeCommits}"

# Create the snapshot commit
cd $_TRILINOS_HOME/cmake/tribits/
../../TriBITS/tribits/snapshot_tribits.py --clean-ignored-files-orig-dir "$@" \
  || echo "snapshot_tribits.py failed but it worked okay?"
cd - 2>&1 >> /dev/null

# Get the current commit message created by the snapshotting script
currentCommitMsg=$(git log -1 --pretty="%B")
#echo "currentCommitMsg:
#
#${currentCommitMsg}"

# Amend the commit message

echo
echo "Amending the snapshot commit message with the new first-parent commits from TriBITS repo:"
echo
echo "${newTribitsMergeCommits}"
echo

newCommitMsg="${currentCommitMsg}

New first-parent commits in TriBITS snapshot:

${newTribitsMergeCommits}\n"

#echo "CommitMsg:
#
#${newCommitMsg}"

git commit -s --amend -m "${newCommitMsg}"
