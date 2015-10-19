#!/bin/bash
#
# Script: dummy_test_commit.sh
#
# This script simulates a test script used with 'git bisect run <script>' to
# show how to use the is_checkin_tested_commit.py script to skip commits that
# are not known to be tested with the checkin-test.py script.  To use this
# script, set the env variable DUMMY_TEST_COMMIT_BAD_SHA to the SHA1 of a
# commit that you are pretending is the bad commit in the range of commits
# <good-commit>..<bad-commit> and then run:
#
#   git bisect <bad-commit> <good-commit>
#   git bisect run ./dummy_test_commit.sh
#
# This should result in git-bisect bounding the commits around
# $DUMMY_TEST_COMMIT_BAD_SHA.  To see the sorted set of commits containig the
# first bad commit, run:
#
#   git bisect log | grep "possible first bad commit"
#

LOG_DUMMY_COMMIT=`git log --oneline HEAD ^$DUMMY_TEST_COMMIT_BAD_SHA^`
if [ "$LOG_DUMMY_COMMIT" == "" ] ; then
  echo "Commit is *before* bad commit $DUMMY_TEST_COMMIT_BAD_SHA!"
else
  echo "Commit is or after *after* bad commit $DUMMY_TEST_COMMIT_BAD_SHA!"
fi    

# Skip the commit if not tested with checkin-test.py script
$TRIBITS_DIR/ci_support/is_checkin_tested_commit.py
IS_CHECKIN_TESTED_COMMIT_RTN=$?
if [ "$IS_CHECKIN_TESTED_COMMIT_RTN" != "0" ] ; then
  exit 125 # Skip the commit because it was not known to be tested!
fi

echo "Building the current version ..."
if [ "$LOG_DUMMY_COMMIT" == "" ] ; then
echo "Commit is *before* bad commit $DUMMY_TEST_COMMIT_BAD_SHA! so marking good!"
  exit 0
else
  echo "Commit is or *after* bad commit $DUMMY_TEST_COMMIT_BAD_SHA so marking bad!"
  exit 1
fi    
