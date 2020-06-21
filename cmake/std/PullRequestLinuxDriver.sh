#!/usr/bin/env bash
# set -x  # echo commands

#    make  sure we have a newer version of git/python
cuda_regex=".*(_cuda_).*"
ride_regex=".*(ride).*"
if [[ ${JOB_BASE_NAME:?} =~ ${cuda_regex} ]]; then
    if [[ ${NODE_NAME:?} =~ ${ride_regex} ]]; then
        echo -e "Job is CUDA"
        module load git/2.10.1
        module load python/2.7.12
    else
        echo -e "ERROR: Unable to find matching environment for CUDA job not on Ride."
        exit -1
    fi
else
    source /projects/sems/modulefiles/utils/sems-modules-init.sh
    module load sems-git/2.10.1
    module load sems-python/2.7.9
fi

# Identify the path to this script
SCRIPTPATH="$(cd "$(dirname "$0")" ; pwd -P)"
echo -e "SCRIPTPATH: ${SCRIPTPATH}"

# Identify the path to the trilinos repository root
REPO_ROOT=`readlink -f ${SCRIPTPATH}/../..`
echo -e "REPO_ROOT : ${REPO_ROOT}"

# This is the old Linux Driver (deprecated)
#${SCRIPTPATH}/PullRequestLinuxDriver-old.sh

# Both scripts will need access through the sandia
# proxies so set them here.
export https_proxy=http://wwwproxy.sandia.gov:80
export http_proxy=http://wwwproxy.sandia.gov:80
export no_proxy='localhost,localnets,127.0.0.1,169.254.0.0/16,forge.sandia.gov'


# Call the script to handle merging the incoming branch into
# the current trilinos/develop branch for testing.
${SCRIPTPATH}/PullRequestLinuxDriverMerge.py ${TRILINOS_SOURCE_REPO:?} \
                                              ${TRILINOS_SOURCE_BRANCH:?} \
                                              ${TRILINOS_TARGET_REPO:?} \
                                              ${TRILINOS_TARGET_BRANCH:?} \
                                              ${TRILINOS_SOURCE_SHA:?} \
                                              ${WORKSPACE:?}

# Call the script to handle driving the testing
${SCRIPTPATH}/PullRequestLinuxDriverTest.py ${TRILINOS_SOURCE_REPO:?} \
                                            ${TRILINOS_SOURCE_BRANCH:?} \
                                            ${TRILINOS_TARGET_REPO:?} \
                                            ${TRILINOS_TARGET_BRANCH:?} \
                                            ${JOB_BASE_NAME:?} \
                                            ${PULLREQUESTNUM:?} \
                                            ${BUILD_NUMBER:?} \
                                            ${WORKSPACE:?}

