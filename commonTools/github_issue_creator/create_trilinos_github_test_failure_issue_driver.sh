#!/bin/bash
#
# Run to create a new Triline GitHub issue body text for nonpassing Trilinos
# GenConfig build test failures shown on CDash.
#
# Usage:
#
#   create_trilinos_github_test_failure_issue_driver.sh \
#     -u "<nonpassing-tests-cdash-url>" \
#     [-s "<summary-line>"] \
#     <extra-args?>
#
# See below implementation for more details

# Find directory path

if [ "$TRILINOS_ATDM_STATUS_DIR" == "" ] ; then
  _ABS_FILE_PATH=`readlink -f $0` || \
    echo "Could not follow symlink to set TRILINOS_ATDM_STATUS_DIR!"
  if [ "$_ABS_FILE_PATH" != "" ] ; then
    TRILINOS_ATDM_STATUS_DIR=`dirname $_ABS_FILE_PATH`
  fi
fi

# Run the tool

time ${TRILINOS_ATDM_STATUS_DIR}/create_trilinos_github_test_failure_issue.py \
  --new-issue-tracker-file=newGithubMarkdownIssueBody.md \
  "$@"
