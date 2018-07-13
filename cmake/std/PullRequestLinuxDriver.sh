#!/bin/env bash

export https_proxy=https://wwwproxy.sandia.gov:80
export http_proxy=http://wwwproxy.sandia.gov:80

whoami
which -a env

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

ls
cd Trilinos

#------------------------------
# Doing merge of pull request
#------------------------------


# Check for existence of source_remote and remove if it exists
git_remote_text=`git remote -v | grep "source_remote"`
if [[ "$git_remote_text" != "" ]]; then
  echo "git remote exists, removing it."
  git remote rm source_remote
fi

#Add the necessary remote
git remote add source_remote ${TRILINOS_SOURCE_REPO:?}
ierror=$?
if [[ $ierror != 0 ]]; then
echo "There was a problem adding the remote for the source repo. The error code was: $ierror"
#just in case somehow a previously defined source_remote caused this failure
#would be better to check prior to the add. Don't want to issue command that will be known to fail typically.
#git remote rm source_remote
exit $ierror
fi

git fetch source_remote
ierror=$?
if [[ $ierror != 0 ]]; then
echo "Source remote fetch failed. The error code was: $ierror"
#git remote rm source_remote
exit $ierror
fi

git remote -v

git merge source_remote/${TRILINOS_SOURCE_BRANCH:?}
ierror=$?
if [[ $ierror != 0 ]]; then
echo "There was an issue merging changes from ${TRILINOS_SOURCE_REPO:?} ${TRILINOS_SOURCE_BRANCH:?} onto ${TRILINOS_TARGET_REPO:?} ${TRILINOS_TARGET_BRANCH:?}. The error code was: $ierror"
#git remote rm source_remote
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
    #git remote rm source_remote
  exit -1
  fi
fi

# may eventually want to verify the same target SHA too, but there would be
# advantages to testing newer versions of target instead of older if
# not all test jobs kick off right away

#------------------------------
# Doing setup for build
#------------------------------


git status
git diff origin/${TRILINOS_TARGET_BRANCH:?} --numstat > ../gitchanges.txt
ierror=$?
if [[ $ierror != 0 ]]; then
echo "There was an issue getting the list of changed files. The error code was: $ierror"

exit $ierror
fi

cd ../
#process list of changes here with script and save to 'packageEnables'
declare packageEnables=$(bash Trilinos/commonTools/test/utilities/changedPackages.bash)
ierror=$?
if [[ $ierror != 0 ]]; then
echo "There was an issue creating the list of package enables. The error code was: $ierror"
exit $ierror
fi
echo "List of package enables is: $packageEnables"

# Set up the full environment for the build
if [ "Trilinos_pullrequest_gcc_4.8.4" == "${JOB_BASE_NAME:?}" ]
then
  source Trilinos/cmake/std/sems/PullRequestGCC4.8.4TestingEnv.sh
  ierror=$?
  if [[ $ierror != 0 ]]; then
    echo "There was an issue loading the gcc environment. The error code was: $ierror"
    exit $ierror
  fi
elif [ "Trilinos_pullrequest_gcc_4.9.3" == "${JOB_BASE_NAME:?}" ]
then
  source Trilinos/cmake/std/sems/PullRequestGCC4.9.3TestingEnv.sh
  ierror=$?
  if [[ $ierror != 0 ]]; then
    echo "There was an issue loading the gcc environment. The error code was: $ierror"
    exit $ierror
  fi
elif [ "Trilinos_pullrequest_intel_17.0.1" == "${JOB_BASE_NAME:?}" ]
then
  source Trilinos/cmake/std/sems/PullRequestIntel17.0.1TestingEnv.sh
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

echo "MPI type = sems-${SEMS_MPI_NAME:?}/${SEMS_MPI_VERSION:?}"

CDASH_TRACK="Pull Request"
echo "CDash Track = ${CDASH_TRACK:?}"


#-------------------------------------
# Doing configure/build/test/submit
#-------------------------------------
echo $packageEnables | sed -e "s/-D\([^= ]*\)=\([^ ]*\)/set(\1 \2 CACHE BOOL \"Enabled by PR package enable file.\")^/g" | tr "^" "\n" > packageEnables.cmake

build_name="PR-$PULLREQUESTNUM-test-$JOB_BASE_NAME-$BUILD_NUMBER"

#This should be runnable from anywhere, but all the tests so far have been from the 
#same dir the simple_testing.cmake file was in.
cd TFW_testing_single_configure_prototype

if [ "icc" == ${CC:?} ]
then
  CONFIG_SCRIPT=PullRequestLinuxIntelTestingSettings.cmake
else
  if [ "Trilinos_pullrequest_gcc_4.8.4" == "${JOB_BASE_NAME:?}" ]
  then
    CONFIG_SCRIPT=PullRequestLinuxGCC4.8.4TestingSettings.cmake
  elif [ "Trilinos_pullrequest_gcc_4.9.3" == "${JOB_BASE_NAME:?}" ]
  then
    CONFIG_SCRIPT=PullRequestLinuxGCC4.9.3TestingSettings.cmake
  fi
fi

ctest -S simple_testing.cmake \
  -Dbuild_name=${build_name:?} \
  -Dskip_by_parts_submit=OFF \
  -Dskip_update_step=ON \
  -Ddashboard_model=Experimental \
  -Ddashboard_track="${CDASH_TRACK:?}" \
  -DPARALLEL_LEVEL=22 \
  -Dbuild_dir="${WORKSPACE:?}/pull_request_test" \
  -Dconfigure_script=../Trilinos/cmake/std/${CONFIG_SCRIPT:?} \
  -Dpackage_enables=../packageEnables.cmake \
  -Dsubprojects_file=../TFW_single_configure_support_scripts/package_subproject_list.cmake

ierror=$?

if [[ $ierror != 0 ]]; then
echo "Single configure/build/test failed. The error code was: $ierror"
exit $ierror
fi


ierror=$?
if [[ $ierror != 0 ]]; then
echo "There was an error removing the source remote. The error code was: $ierror"
exit $ierror
fi

#NEED TO MAKE SURE THE REPOS ARE CLEAN FOR NEW PULL REQUESTS!

#pushd Trilinos/cmake/ctest/drivers/parameterized
#ctest -S ctest_linux_nightly_generic.cmake
