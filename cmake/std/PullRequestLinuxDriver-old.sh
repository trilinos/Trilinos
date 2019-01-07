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

export https_proxy=http://wwwproxy.sandia.gov:80
export http_proxy=http://wwwproxy.sandia.gov:80
no_proxy='localhost,localnets,.sandia.gov,127.0.0.1,169.254.0.0/16,forge.sandia.gov'

whoami
which -a env

echo -e ""
echo -e "Jenkins Environment Variables:"
echo -e "- JOB_BASE_NAME: ${JOB_BASE_NAME}"
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
module load sems-git/2.10.1


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
cd Trilinos

# Check for existence of source_remote and remove if it exists
git_remote_text=`git remote -v | grep "source_remote"`
if [[ "$git_remote_text" != "" ]]; then
  echo "git remote exists, removing it."
  git remote rm source_remote
fi

# Add the necessary remote
git remote add source_remote ${TRILINOS_SOURCE_REPO:?}
ierror=$?
if [[ $ierror != 0 ]]; then
  echo "There was a problem adding the remote for the source repo. The error code was: $ierror"
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
    echo "Source remote fetch failed. The error code was: $ierror"
    if [[ $i != $num_retries ]]
    then
      echo "retry $i"
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
  echo "Origin target remote fetch failed. The error code was: $ierror"
  exit $ierror
fi

git status

# Clean out any non-committed changes so that we can checkout the correct
# branch
git reset --hard HEAD
ierror=$?
if [[ $ierror != 0 ]]; then
  echo "There was an error clearing out any local uncommitted changes. The error code was: $ierror"
  exit $ierror
fi

git status

# Get on the right local branch
git checkout -B ${TRILINOS_TARGET_BRANCH:?} origin/${TRILINOS_TARGET_BRANCH:?}
ierror=$?
if [[ $ierror != 0 ]]; then
  echo "There was an error checking out and updating to the target remote branch. The error code was: $ierror"
  exit $ierror
fi

git status

# Merge the souce branch into the local target branch
git merge --no-edit source_remote/${TRILINOS_SOURCE_BRANCH:?}
ierror=$?
if [[ $ierror != 0 ]]; then
  echo "There was an issue merging changes from ${TRILINOS_SOURCE_REPO:?} ${TRILINOS_SOURCE_BRANCH:?} onto ${TRILINOS_TARGET_REPO:?} ${TRILINOS_TARGET_BRANCH:?}. The error code was: $ierror"
  exit $ierror
fi

#Need to compare expected SOURCE SHA to actual SHA! This will prevent a security hole.
#first get the most recent SHA on the source branch
declare source_sha=$(git rev-parse source_remote/${TRILINOS_SOURCE_BRANCH:?})
echo "The most recent SHA for repo: ${TRILINOS_SOURCE_REPO:?} on branch: ${TRILINOS_SOURCE_BRANCH:?} is: $source_sha"
#Now see if the two shas match, unless TRILINOS_SOURCE_SHA is the default value of ABC
if [[ ABC != ${TRILINOS_SOURCE_SHA:?} ]]; then
  if [[ $source_sha != ${TRILINOS_SOURCE_SHA:?} ]]; then
    echo "The SHA ($source_sha) for the last commit on branch ${TRILINOS_SOURCE_BRANCH:?} in repo ${TRILINOS_SOURCE_REPO:?} is different than the expected SHA, which is: ${TRILINOS_SOURCE_SHA:?}. The error code was: $ierror"
  exit -1
  fi
fi

# may eventually want to verify the same target SHA too, but there would be
# advantages to testing newer versions of target instead of older if
# not all test jobs kick off right away

#------------------------------
# PR merge is complete
#------------------------------

# TODO: Split this script in two here (See Issue 3625 for reasons)

#------------------------------
# Doing setup for build
#------------------------------

cd ../

# Set up the full environment for the build
if [ "Trilinos_pullrequest_gcc_4.8.4" == "${JOB_BASE_NAME:?}" ] ; then
  source ${TRILINOS_DRIVER_SRC_DIR}/cmake/std/sems/PullRequestGCC4.8.4TestingEnv.sh
  ierror=$?
  if [[ $ierror != 0 ]]; then
    echo "There was an issue loading the gcc environment. The error code was: $ierror"
    exit $ierror
  fi
elif [ "Trilinos_pullrequest_gcc_4.9.3" == "${JOB_BASE_NAME:?}" ] ; then
  source ${TRILINOS_DRIVER_SRC_DIR}/cmake/std/sems/PullRequestGCC4.9.3TestingEnv.sh
  ierror=$?
  if [[ $ierror != 0 ]]; then
    echo "There was an issue loading the gcc environment. The error code was: $ierror"
    exit $ierror
  fi
elif [ "Trilinos_pullrequest_gcc_4.9.3_SERIAL" == "${JOB_BASE_NAME:?}" ] ; then
  # TODO: Update this to use a 4.9.3 SERIAL testing environment script.
  source ${TRILINOS_DRIVER_SRC_DIR}/cmake/std/sems/PullRequestGCC4.9.3TestingEnvSERIAL.sh 
  ierror=$?
  if [[ $ierror != 0 ]]; then
    echo "There was an issue loading the gcc environment. The error code was: $ierror"
    exit $ierror
  fi
elif [ "Trilinos_pullrequest_intel_17.0.1" == "${JOB_BASE_NAME:?}" ] ; then
  source ${TRILINOS_DRIVER_SRC_DIR}/cmake/std/sems/PullRequestIntel17.0.1TestingEnv.sh
  ierror=$?
  if [[ $ierror != 0 ]]; then
    echo "There was an issue loading the intel environment. The error code was: $ierror"
    exit $ierror
  fi
else
  ierror=42
  echo "There was an issue loading the proper environment. The error code was: $ierror"
  exit $ierror
fi

# The single submit requires at least cmake 3.10.*
cmake --version

module list

# This crashes for the serial case since MPI variables are not set
# - See Issue #3625
# - wcm: bugfix #3673
regex=".*(_SERIAL)$"
if [[ ! ${JOB_BASE_NAME:?} =~ ${regex} ]]; then
  echo "MPI type = sems-${SEMS_MPI_NAME:?}/${SEMS_MPI_VERSION:?}"
else
  echo "Job is SERIAL"
fi

CDASH_TRACK="Pull Request"
echo "CDash Track = ${CDASH_TRACK:?}"


#-------------------------------------
# Doing configure/build/test/submit
#-------------------------------------

rm packageEnables.cmake
env \
  TRILINOS_DIR=${PR_TRILINOS_DIR} \
  TRILINOS_SCRIPTS_DIR=${TRILINOS_DRIVER_SRC_DIR} \
${TRILINOS_DRIVER_SRC_DIR}/commonTools/framework/get-changed-trilinos-packages.sh \
  origin/${TRILINOS_TARGET_BRANCH:?} HEAD packageEnables.cmake
# NOTE: Above we use the git diff origin/<target-branch>..HEAD to give us the
# correct list of changed files.  This works because this is done after
# merging the target branch and the soruce branch.  With that, git diff will
# show all of the changes in the merged copy from what is in <target-branch>.
# This with give the correct set of changed files even if older versions of
# <target-branch> were merged multiple times into <source-branch>.  It turns
# out that the only way to have git show the correct set of diffs in that case
# is to actually do the merge and then do the diff.

ierror=$?
if [[ $ierror != 0 ]]; then
  echo "There was an issue generating packageEnables.cmake.  The error code was: $ierror"
  exit $ierror
fi

echo "Enabled packages:"
cmake -P packageEnables.cmake

build_name="PR-$PULLREQUESTNUM-test-$JOB_BASE_NAME-$BUILD_NUMBER"

#This should be runnable from anywhere, but all the tests so far have been from the 
#same dir the simple_testing.cmake file was in.
cd TFW_testing_single_configure_prototype

if [ "icc" == ${CC:?} ] ; then
  CONFIG_SCRIPT=PullRequestLinuxIntelTestingSettings.cmake
else
  if [ "Trilinos_pullrequest_gcc_4.8.4" == "${JOB_BASE_NAME:?}" ]; then
    CONFIG_SCRIPT=PullRequestLinuxGCC4.8.4TestingSettings.cmake
  elif [ "Trilinos_pullrequest_gcc_4.9.3" == "${JOB_BASE_NAME:?}" ]; then
    CONFIG_SCRIPT=PullRequestLinuxGCC4.9.3TestingSettings.cmake
  elif [ "Trilinos_pullrequest_gcc_4.9.3_SERIAL" == "${JOB_BASE_NAME:?}" ]; then
    # TODO: Update this to use a 4.9.3 SERIAL testing environment script.
    CONFIG_SCRIPT=PullRequestLinuxGCC4.9.3TestingSettingsSERIAL.cmake
  fi
fi

ctest -S simple_testing.cmake \
  -Dbuild_name=${build_name:?} \
  -Dskip_by_parts_submit=OFF \
  -Dskip_update_step=ON \
  -Ddashboard_model=Experimental \
  -Ddashboard_track="${CDASH_TRACK:?}" \
  -DPARALLEL_LEVEL=18 \
  -Dbuild_dir="${WORKSPACE:?}/pull_request_test" \
  -Dconfigure_script=${TRILINOS_DRIVER_SRC_DIR}/cmake/std/${CONFIG_SCRIPT:?} \
  -Dpackage_enables=../packageEnables.cmake \
  -Dsubprojects_file=../TFW_single_configure_support_scripts/package_subproject_list.cmake

ierror=$?
if [[ $ierror != 0 ]]; then
  echo "Single configure/build/test failed. The error code was: $ierror"
  exit $ierror
fi

#pushd Trilinos/cmake/ctest/drivers/parameterized
#ctest -S ctest_linux_nightly_generic.cmake
