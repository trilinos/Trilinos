#!/usr/bin/env bash
# set -x  # echo commands

# Identify the path to this script
SCRIPTPATH="$(cd "$(dirname "$0")" ; pwd -P)"
echo -e "SCRIPTPATH: ${SCRIPTPATH}"


# Both scripts will need access through the sandia
# proxies so set them here.
export https_proxy=http://wwwproxy.sandia.gov:80
export http_proxy=http://wwwproxy.sandia.gov:80
export no_proxy='localhost,localnets,127.0.0.1,169.254.0.0/16,forge.sandia.gov'


# Call the script to handle merging the incoming branch into
# the current trilinos/develop branch for testing.
${SCRIPTPATH}/PullRequestLinuxDriverMerge.py ${TRILINOS_SOURCE_REPO:?}    \
                                             ${TRILINOS_SOURCE_BRANCH:?}  \
                                             ${TRILINOS_TARGET_REPO:?}    \
                                             ${TRILINOS_TARGET_BRANCH:?}  \
                                             ${TRILINOS_SOURCE_SHA:?}     \
                                             ${WORKSPACE:?}

err=$?
if [ $err != 0 ]; then
    exit $err
fi

