#!/usr/bin/env bash
# set -x  # echo commands

# Identify the path to this script
SCRIPTPATH="$(cd "$(dirname "$0")" ; pwd -P)"
echo -e "SCRIPTPATH: ${SCRIPTPATH}"

# Identify the path to the trilinos repository root
REPO_ROOT=`readlink -f ${SCRIPTPATH}/../..`
echo -e "REPO_ROOT : ${REPO_ROOT}"

# This is the old Linux Driver (deprecated)
#${SCRIPTPATH}/PullRequestLinuxDriver-old.sh

# Call the script to handle merging the incoming branch into
# the current trilinos/develop branch for testing.
${SCRIPTPATH}/PullRequestLinuxDriver-Merge.sh

# Call the script to handle driving the testing
${SCRIPTPATH}/PullRequestLinuxDriver-Test.sh

