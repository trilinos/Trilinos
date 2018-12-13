#!/usr/bin/env bash
# set -x  # echo commands

#
# This script drives a PR testing build.  It assume that Trilinos is already
# cloned under $PWD/Trilinos and that the 'origin' remote points to
# $TRILINOS_TARGET_REPO (but that is not checked here).
#
# As long as the ${PWD}/Trilinos git repo has the correct 'origin', this
# script will automatically set it up to do the merge correctly, no matter
# what its state before this script is called (i.e. from a past PR
# attempt). Unless the Trilinos/.git directory becomes corrupted, there should
# *NEVER* be any need to delete and reclone this Trilinos git repo.
#
# This script can be run in a mode where the driver scripts are run from one
# Trilinos git repo and operate on another Trilinos git repo that gets
# manipulated in order to merge the "source" topic branch into the "target"
# branch.  This makes it easy to test changes to the PR scripts.  But if this
# script is run from ${PWD}/Trilinos, then these repos are one and the same
# and we get the correct behavior for PR testing.
#
echo -e "--------------------------------------------------------------------------------"
echo -e "-"
echo -e "- Begin: PullRequestLinuxDriver-Merge.sh"
echo -e "-"
echo -e "--------------------------------------------------------------------------------"

# This script expects to start out in the root level of the Jenkins workspace.
# Let's make sure we're there.
cd ${WORKSPACE:?}

# Identify the path to this script
SCRIPTPATH="$(cd "$(dirname "$0")" ; pwd -P)"
echo -e "SCRIPTPATH: ${SCRIPTPATH:?}"

# Identify the path to the trilinos repository root
REPO_ROOT=`readlink -f ${SCRIPTPATH}/../..`
echo -e "REPO_ROOT : ${REPO_ROOT:?}"

# Set Sandia Proxy environment variables
export https_proxy=http://wwwproxy.sandia.gov:80
export http_proxy=http://wwwproxy.sandia.gov:80
no_proxy='localhost,localnets,.sandia.gov,127.0.0.1,169.254.0.0/16,forge.sandia.gov'

# Print out some useful environment information
whoami
which -a env

# Useful information for when it counts.
echo -e ""
echo -e "================================================================================"
echo -e "Jenkins Environment Variables:"
echo -e "- JOB_BASE_NAME: ${JOB_BASE_NAME:?}"
echo -e "- JOB_NAME     : ${JOB_NAME:?}"
echo -e "- WORKSPACE    : ${WORKSPACE:?}"
echo -e "- NODE_NAME    : ${NODE_NAME:?}"
echo -e ""
echo -e "================================================================================"
echo -e "Environment:"
echo -e ""
echo -e "  pwd = `pwd`"
echo -e ""
env
echo -e ""
echo -e "================================================================================"
echo -e ""

## Rather than do proper option handling right now I am just going to
##  test that all these environment variables are set.  Because getopt ;->
: ${TRILINOS_SOURCE_REPO:?}
: ${TRILINOS_SOURCE_BRANCH:?}
: ${TRILINOS_TARGET_REPO:?}
: ${TRILINOS_TARGET_BRANCH:?}
: ${TRILINOS_SOURCE_SHA:?}
: ${PULLREQUESTNUM:?}
: ${JOB_BASE_NAME:?}
: ${BUILD_NUMBER:?}
: ${WORKSPACE:?}

source /projects/sems/modulefiles/utils/sems-modules-init.sh

declare -i ierror=0
#Have to keep loading git
regex=".*(_cuda_).*"
if [[ ! ${JOB_BASE_NAME:?} =~ ${regex} ]]; then
  module load sems-git/2.10.1
else
  echo -e "Job is CUDA and assumed on Ride"
  module load git/2.10.1
fi

#--------------------------------------------
# Get Trilinos scripts and PR merge repo dirs
#--------------------------------------------

# The Trilinos dir that the PR merge will be done in
PR_TRILINOS_DIR=$PWD/Trilinos
echo "PR_TRILINOS_DIR = $PR_TRILINOS_DIR"

# Get the Trilinos scripts driver dir where this script is run from.  First,t
# ry to grab from the symlink (only works on Linux)
_ABS_FILE_PATH=`readlink -f $0` || \
    echo "Could not follow symlink to set TRILINOS_DRIVER_SRC_DIR!"

if [ "$_ABS_FILE_PATH" != "" ] ; then
    _SCRIPT_DIR=`dirname $_ABS_FILE_PATH`
    TRILINOS_DRIVER_SRC_DIR=$_SCRIPT_DIR/../..
else
    # If that did not work, then we are not on Linux so give up and assume the
    # standard location
    TRILINOS_DRIVER_SRC_DIR=${PR_TRILINOS_DIR}
fi
echo "TRILINOS_DRIVER_SRC_DIR = $TRILINOS_DRIVER_SRC_DIR"


#------------------------------
# Doing merge of pull request
#------------------------------

ls
pushd Trilinos &> /dev/null
echo -e "Set CWD = `pwd`"

# Check for existence of source_remote and remove if it exists
git_remote_text=`git remote -v | grep "source_remote"`
if [[ "$git_remote_text" != "" ]]; then
    echo -e "git remote exists, removing it."
    git remote rm source_remote
fi

# Add the necessary remote
git remote add source_remote ${TRILINOS_SOURCE_REPO:?}
ierror=$?
if [[ $ierror != 0 ]]; then
    echo -e "There was a problem adding the remote for the source repo. The error code was: $ierror"
    #just in case somehow a previously defined source_remote caused this failure
    #would be better to check prior to the add. Don't want to issue command that
    #will be known to fail typically.  git remote rm source_remote
    exit $ierror
fi

git remote -v

num_retries=3

for i in `seq ${num_retries}`
do
    git fetch source_remote ${TRILINOS_SOURCE_BRANCH:?}
    ierror=$?
    if [[ $ierror != 0 ]]; then
        echo -e "Source remote fetch failed. The error code was: $ierror"
        if [[ $i != $num_retries ]]; then
            echo -e "retry $i"
            sleep $(($i*20))
        else
            exit $ierror
        fi
    else
        break
    fi
done

git fetch origin ${TRILINOS_TARGET_BRANCH:?}
ierror=$?
if [[ $ierror != 0 ]]; then
    echo -e "Origin target remote fetch failed. The error code was: $ierror"
    exit $ierror
fi

git status

# Clean out any non-committed changes so that we can checkout the correct
# branch
git reset --hard HEAD
ierror=$?
if [[ $ierror != 0 ]]; then
    echo -e "There was an error clearing out any local uncommitted changes. The error code was: $ierror"
    exit $ierror
fi

git status

# Get on the right local branch
git checkout -B ${TRILINOS_TARGET_BRANCH:?} origin/${TRILINOS_TARGET_BRANCH:?}
ierror=$?
if [[ $ierror != 0 ]]; then
    echo -e "There was an error checking out and updating to the target remote branch. The error code was: $ierror"
    exit $ierror
fi

git status

# Merge the souce branch into the local target branch
git merge --no-edit source_remote/${TRILINOS_SOURCE_BRANCH:?}
ierror=$?
if [[ $ierror != 0 ]]; then
    echo -e "There was an issue merging changes from "
    echo -e "  ${TRILINOS_SOURCE_REPO:?} ${TRILINOS_SOURCE_BRANCH:?} onto ${TRILINOS_TARGET_REPO:?} ${TRILINOS_TARGET_BRANCH:?}."
    echo -e "  The error code was: $ierror"
    exit $ierror
fi

#Need to compare expected SOURCE SHA to actual SHA! This will prevent a security hole.
#first get the most recent SHA on the source branch
declare source_sha=$(git rev-parse source_remote/${TRILINOS_SOURCE_BRANCH:?})
echo "The most recent SHA for repo: ${TRILINOS_SOURCE_REPO:?} on branch: ${TRILINOS_SOURCE_BRANCH:?} is: $source_sha"
#Now see if the two shas match, unless TRILINOS_SOURCE_SHA is the default value of ABC
if [[ ABC != ${TRILINOS_SOURCE_SHA:?} ]]; then
    if [[ $source_sha != ${TRILINOS_SOURCE_SHA:?} ]]; then
        echo -e "The SHA ($source_sha) for the last commit on branch ${TRILINOS_SOURCE_BRANCH:?}"
        echo -e "  in repo ${TRILINOS_SOURCE_REPO:?} is different than the expected SHA,"
        echo -e "  which is: ${TRILINOS_SOURCE_SHA:?}. The error code was: $ierror"
        exit -1
    fi
fi

# may eventually want to verify the same target SHA too, but there would be
# advantages to testing newer versions of target instead of older if
# not all test jobs kick off right away

#------------------------------
# PR merge is complete
#------------------------------

# Return to previous directory
popd &> /dev/null
echo -e "Set CWD = `pwd`"
